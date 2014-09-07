/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Note: Algorithms are implemented from paper: "Simplification and extension of Spread Constraint by Pierre et al" 

/**@file   cons_Spread.c
 * @brief  constraint handler for Spread constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

//#define SCIP_DEBUG
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "cons_spread.h"
#include "scip/cons_linking.h"
#include "scip/cons_knapsack.h"
#include "scip/scipdefplugins.h"

//headers used by bound consistency algorithm
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <algorithm>
#include <math.h>


//fundamental constraint handler properties 
#define CONSHDLR_NAME          "Spread"
#define CONSHDLR_DESC          "Spread constraint handler"
#define CONSHDLR_SEPAPRIORITY   2100000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -3000000 //< priority of the constraint handler for constraint enforcing 
#define CONSHDLR_CHECKPRIORITY -3000000 //< priority of the constraint handler for checking feasibility 
#define CONSHDLR_EAGERFREQ          100 //< frequency for using all instead of only the useful constraints in separation,                                                 propagation and enforcement, -1 for no eager evaluations, 0 for first only 
// Do you need this constraint handler? if no constraint handler available ?
#define CONSHDLR_NEEDSCONS         TRUE //< should the constraint handler be skipped, if no constraints are available? 

// optional constraint handler properties 
// TODO: remove properties which are never used because the corresponding routines are not supported 
#define CONSHDLR_SEPAFREQ            1 //< frequency for separating cuts; zero means to separate only in the root node 
#define CONSHDLR_DELAYSEPA        FALSE //< should separation method be delayed, if other separators found cuts? 
#define CONSHDLR_PROPFREQ            1 //< frequency for propagating domains; zero means only preprocessing propagation 
#define CONSHDLR_DELAYPROP        FALSE //< should propagation method be delayed, if other propagators found reductions? 
#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP//< propagation timing mask of the constraint handler
#define CONSHDLR_MAXPREROUNDS        -1 //< maximal number of presolving rounds the constraint handler participates in (-1: no limit) 
#define CONSHDLR_DELAYPRESOL      FALSE //< should presolving method be delayed, if other presolvers found reductions? 

// (optional) enable linear or nonlinear constraint upgrading 
#include "scip/cons_linear.h"
#include "scip/cons_nonlinear.h"

/* default parameter values */

/* separation */
#define DEFAULT_USEBINVARS             FALSE //< should the binary representation be used? 
#define DEFAULT_LOCALCUTS              FALSE //< should cuts be added only locally? 
#define DEFAULT_USECOVERCUTS            TRUE //< should covering cuts be added? 
#define DEFAULT_CUTSASCONSS             TRUE //< should the cuts be created as knapsack constraints? 
#define DEFAULT_SEPAOLD                 TRUE //< shall old sepa algo be applied? 

// propagation 
#define DEFAULT_CORETIMES               TRUE //< should core-times be propagated (time tabling)? 
#define DEFAULT_OVERLOAD                TRUE //< should edge finding be used to detect an overload? 
#define DEFAULT_EDGEFINDING             TRUE //< should edge-finding be executed? 
#define DEFAULT_USEADJUSTEDJOBS        FALSE //< should during edge-finding jobs be adusted which run on the border of the effective time horizon? 

// presolving 
#define DEFAULT_DUALPRESOLVE            TRUE //< should dual presolving be applied? 
#define DEFAULT_COEFTIGHTENING         FALSE //< should coefficient tightening be applied? 
#define DEFAULT_NORMALIZE               TRUE //< should demands and capacity be normalized? 
#define DEFAULT_MAXNODES             10000LL //< number of branch-and-bound nodes to solve an independent Spread constraint  (-1: no limit) 

// enforcement 
#define DEFAULT_FILLBRANCHCANDS        FALSE //< should branching candidates be added to storage? 

// conflict analysis 
#define DEFAULT_USEBDWIDENING           TRUE //< should bound widening be used during conflict analysis? 

#define EVENTHDLR_NAME         "Spread"
#define EVENTHDLR_DESC         "bound change event handler for Spread constraints"
/*
Hi, I did not have time to think about these and I realize this today on my last day on research.
So, when testing, if you encounter errors, the errors might be at the below cases.
Important ReadMe:
1. For the general mean value, I pick the q that resulted in the largest
DMAX(q) function. Note: I did not account for the cases where q results in shifting to a previous interval,
and therefore need to calculate q recursively for every q at each intervals end points.
The change can be done easily, but I am not sure if the current code is correct,
or I am suppose to calculate the largest q from DMAX(q) function recursively to select q.

2. There is an error in the code. I prune each X individually before I updated them.
That shouldn't be right as it would end up pruning more than allowed based on given stdDevUB.
The right method should be break out of the current iteration as soon as 1 X changes bounds,
and to re-iterate from the beginning. This happens in my code for both the specific and general mean values case.
You can find it by searching for TODO
All the TODO comments are things that needs to be done that I haven't completed.

3. There is a part of the code under TODO: that I did not account for m being equal to 0. Currently,
it does not cause any problems from my test cases, but if you start seeing infinite, nan or negative
values in the bounds, this should be the problem.

*/

// This file contains my implementation of Spread Constraint
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <algorithm>
#include <math.h>

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdarg.h>  // to accept variable arguments
using namespace std;

// Note: The paper does not account for the special case where
// i) m = 0, if m = 0, lots of equations divides by m, division by 0 is undefined
// ii) m = n, if m = n, some equatons divides by (n-m), division by 0 is undefined

// m == n is used to calculate q0 for opt(q), which is needed to prune U using S.
// However, can be proven if m == n, the q0 won't be in that interval, so no worries.

// m == 0  is used to calculate DMax(q)
// However, can note that if m = 0, we don't need to calculate v as no X values will be assigned to v
// so we re-calculate everything by ignoring v, in other words, set v = 0.
// this results in

// if (m == 0)
//{
//    a = 1.0;
//    v = 0;
//    opt doesnt add the ((q-ES)^2/m) term
//}

// else (normal case)
// a = 1.0 + 1.0/m
// v = (q-ES)/m
// opt adds the ((q-ES)^2/m) term

// Ignore all
class Interval
{
    friend class domainVar; // allows domainVar to easily access all public and private members of this class
public: // to make life easier, just call these public,
        // so that main function that contains propagation algorithm can access these easily.
    double Ilb; // interval lower bound
    double Iub; // interval upper bound
    double Vlb; // q's lower bound in this interval to locate v
    double Vub; // q's upper bound in this interval to locate v
    double ES; // extreme sum for all X in this interval, i.e. min + max
    double C; // extreme^2 sum for all X in this interval, i.e. (min)^2 + (max)^2
    int *M; // array of indexes of all X in the M(Iq)
    int *R; // array of indexes of all X in the R(Iq)
    int *L; // array of indexes of all X in the L(Iq)
    int m; // cardinality of M(Iq)
    int r; // cardinality of R(Iq)
    int l; // cardinality of l(Iq)

    // The range of nVar that can be in this interval
     // Note: nVarUb may be lower than nVarLb
    // It may also not contain the minimum value if it contains the global minimum since nVar is a convex function
    double nVarUb;
    double nVarLb;

    double q0; // the value of q0 if this interval was the global min
                // -1 => m == n for this interval

public:
    Interval()
    {
        Ilb = 0;
        Iub = 0;
    }

    // Constructer, initializes with the upper and lower bound of this class
    Interval(int lb, int ub)
    {
        this->Ilb = lb;
        this->Iub = ub;
    }
};

// This class represents a domainVariable
// It includes additional private data for:
// i) Spread Constraint
class domainVar
{
private:
// General for Every Constraint
    int n; // number of X variables
    int **X; // a 2D array of X's and its upper and lower bound
    char name = 'a';

// mergesort for create interval, returns number of values currently in array
// note: This merge sort is unique as it removes duplicate elements while sorting
template <class X>
int mergesort(X a[], int n)
{
    if (n==1)
    {
        return 1;
    }
    int q, p;
    q = n/2;
    p = n/2 + 1;
    X b[q];
    X c[p];
    int i = 0;
    for (i = 0; i < q; i++)
    {
        b[i] = a[i];
    }
    int k = 0;
    for (int j = q; j < n; j++)
    {
        c[k] = a[j];
        k++;
    }
    q = mergesort(b, i);
    p = mergesort(c, k);
    int r, s, t;
    t = 0; r = 0; s = 0;
    while( (r!=q) && (s!= p))
    {
        if (b[r] < c[s])
        {
            a[t] = b[r];
            r++;
        }
        else if (b[r] > c[s])
        {
            a[t] = c[s];
            s++;
        }
        else // here it means both of them are equal
        {
            a[t] = b[r];
            r++;
            s++; // also incrememnt s to skip it
        }
        t++;

    }
    if ( r==q)
    {
        while(s!=p)
        {
            a[t] = c[s];
            s++;
            t++;
        }
    }
    else
    {
        while(r != q)
        {
            a[t] = b[r];
            r++;
            t++;
        }
    }
    return t;
}

public:
// Specially For Spread Constraint
// Let these variables be public so main function can access these easily.
    int numI; // number of contiguous intervals
    Interval **I; // array of contiguous interval bounds (note: This is a 1D array, but simple hack works with 2 used as 1)
    bool needDelete; // to determine if this is the first time interval I is being created, or needs to be clean up
                    // before creating again

    // constructor
    domainVar()
    {
        n = 10;
        X = new int* [n];
        int i = 0;
        for (i = 0; i < n; i++)
        {
            X[i] = new int [2];
            X[i][0] = 0;
            X[i][1] = 0;
        }
        needDelete = false;
    }

    domainVar(int n_args, char symbol)
    {
        this->n = n_args;
        name = symbol;
        X = new int* [n_args];
        int i = 0;
        for (i = 0; i < n_args; i++)
        {
            X[i] = new int [2];
            X[i][0] = 0;
            X[i][1] = 0;
        }
        needDelete = false;
    }

    void updateDomainVar(int index, int lb, int ub)
    {
        this->setDomVarLb(index, lb);
        this->setDomVarUb(index, ub);
    }

    domainVar(char symbol, int n_args, ...)
    {
        name = symbol;
        va_list head;
        va_start(head, n_args); // initialize head to be every variable after n_args
        n = n_args;
        X = new int* [n];

        int lb, ub;
        for(int i = 0; i < n; i++)
        {
            X[i] = new int [2];
            lb = va_arg(head, int);
            ub = va_arg(head, int);
            X[i][0] = lb;
            X[i][1] = ub;
        }
        va_end(head);
        needDelete = false;
    }

    // Public Methods
    void printDomVar()
    {
        int i = 0;
        for(i = 0; i < n; i++)
        {
            cout << "Bounds of "<< name << "[ " << i << "] : " <<  X[i][0] <<" , " << X[i][1] <<endl;
        }
    }

    void setDomVarLb(int i, int value)
    {
        X[i][0] = value;
    }

    int getDomVarLb(int i)
    {
        return X[i][0];
    }

    void setDomVarUb(int i, int value)
    {
        X[i][1] = value;
    }

    int getDomVarUb(int i)
    {
        return X[i][1];
    }

    void setDomN(int value)
    {
        n = value;
    }

    int getDomN()
    {
        return n;
    }
    double maxStdDevThiago()
    {
        // Maximum std. deviation possible from all values of X
        double maxScalc = 0;
        double sumXi = 0;
        double sumXj = 0;

        double Sub = 0; // Sum of upper bounds
        double Slb = 0; // Sum of lower bounds
        int i = 0;
        for (i = 0; i < this->getDomN(); i++ )
        {
            Slb += this->getDomVarLb(i);
            Sub += this->getDomVarUb(i);
        }

        double minUcalc = Slb/this->getDomN();
        double maxUcalc = Sub/this->getDomN();

   double min_mean_;		//
   double max_mean_;		//
   double distance_x_max;	//
   double distance_x_min;	//
    double X[this->getDomN()];
   // Assign each variable, x[k] to its extrema, such that the std dev is maximized
   // Assigning it to the extrema will maximize the std. dev as explained in the paper
   // "pg 13 of 15  Simplification and extension of Spread Constraint by Pierre et al"
   for (int k=0;k < this->n;k++)
   {
	  // assume x = xmax => min_mean increases
      min_mean_ = minUcalc + (this->getDomVarUb(k) - this->getDomVarLb(k))/this->n ;

// As explained in the paper, find out which values to use before assigning values of x in order to get max std. deviation.
		// lower bound of x = xmax
        if ( fabs(this->getDomVarUb(k) - min_mean_) <= fabs(this->getDomVarUb(k) - maxUcalc) )
        {
            distance_x_max = fabs(this->getDomVarUb(k) - min_mean_);
        }
        else
        {
            distance_x_max = fabs(this->getDomVarUb(k) - maxUcalc);
        }

      //assume x = xmin => max_mean decreases
      max_mean_ = maxUcalc - (this->getDomVarUb(k) - this->getDomVarLb(k))/this->n;
	  // upper bound of x = xmin
	  if (fabs(this->getDomVarLb(k) - max_mean_) >= fabs(this->getDomVarLb(k) - minUcalc))
      {
          distance_x_min = fabs(this->getDomVarLb(k) - max_mean_);
      }
      else
      {
          distance_x_min =  fabs(this->getDomVarLb(k) - minUcalc);
      }

      if (distance_x_max > distance_x_min)
      {
			X[k] = this->getDomVarUb(k); // assign to upper bound extreme
			minUcalc = min_mean_; // update the minimum mean
      }

      else // LINE HAHA THIS SHOULDN'T BE CORRECT
      {
         X[k] = this->getDomVarLb(k); // assign to lower bound extreme
		 maxUcalc = max_mean_; // update the maximum mean
      }
   }
   double std_dev = 0;
   // Here, we're done assigning X to its appropriate extreme values,
    // calculate the maximum std. deviatio that can occur
   for(int i = 0; i < this->n; ++i )
      std_dev += pow(X[i] - 3, 2.0);  // numerator of variance calculation
   std_dev = std_dev/this->n; // variance
   std_dev = sqrt(std_dev); // std_deviation
   return std_dev;
    }

    // This algorithm from Wen Yang is faster but not as tight as the general algorithm given on paper
    double maxStdDevSpecific(double mean)
    {
        double Xi[this->getDomN()];
        double maxScalc = 0;
        double diffOne = 0;
        double diffTwo = 0;
        double sumXi = 0;
        // Assign all Xi
        for (int i = 0; i < this->getDomN(); i++)
        {
            diffOne = fabs(this->getDomVarLb(i) - mean);
            diffTwo = fabs(this->getDomVarUb(i) - mean);
      //      cout <<" diffONe is : " << diffOne<< " and DiffTWO is : " << diffTwo << endl;
            if (diffOne >= diffTwo)
            {
                Xi[i] = this->getDomVarLb(i);
            }
            else
            {
                Xi[i] = this->getDomVarUb(i);
            }
     //       cout << " X[ " << i << " ] is " << Xi[i] << endl;
            sumXi += fabs(Xi[i] - mean);
        }
        maxScalc = sumXi * sumXi;
        maxScalc = maxScalc/this->getDomN();
        maxScalc = sqrt(maxScalc);
        return maxScalc;
    }
    // Calculate and returns max Std Dev from current values of X
    double maxStdDev()
    {
        // Maximum std. deviation possible from all values of X
        double maxScalc = 0;
        double sumXi = 0;
        double sumXj = 0;

        double Sub = 0; // Sum of upper bounds
        double Slb = 0; // Sum of lower bounds
        int i = 0;
        for (i = 0; i < this->getDomN(); i++ )
        {
            Slb += this->getDomVarLb(i);
            Sub += this->getDomVarUb(i);
        }

        double minUcalc = Slb/this->getDomN();
        double maxUcalc = Sub/this->getDomN();

        // To keep track of the terms to substitute in when calculating maxScalc
        double Xi[this->getDomN()];
        double Xj[this->getDomN()];
        double minUcalc_ , maxUcalc_; // _ => The updated values for changing the mean and max to other extreme
        double distXmax[2]; // Lower and Upper bound for distance when X = Xmax
        double distXmin[2]; // Lower and Upper bound for distance when X = Xmin

        for (i = 0; i < this->getDomN(); i++)
        {
            minUcalc_ = minUcalc + ((this->getDomVarUb(i) - this->getDomVarLb(i))/(this->getDomN()));
            maxUcalc_ = maxUcalc - ((this->getDomVarUb(i) - this->getDomVarLb(i))/(this->getDomN()));
            distXmin[0] = std::min(fabs(this->getDomVarLb(i) - minUcalc), fabs(this->getDomVarLb(i) - maxUcalc_));
            distXmin[1] = std::max(fabs(this->getDomVarLb(i) - minUcalc), fabs(this->getDomVarLb(i) - maxUcalc_));
            distXmax[0] = std::min(fabs(this->getDomVarUb(i) - maxUcalc), fabs(this->getDomVarUb(i) - minUcalc_));
            distXmax[1] = std::max(fabs(this->getDomVarUb(i) - maxUcalc), fabs(this->getDomVarUb(i) - minUcalc_));

            if (distXmin[0] >= distXmax[1])
            {
                Xi[i] = this->getDomVarLb(i);
                Xj[i] = this->getDomVarLb(i);
            }

            else if (distXmax[0] >= distXmin[1])
            {
                Xi[i] = this->getDomVarUb(i);
                Xj[i] = this->getDomVarUb(i);
            }

            else if (this->getDomVarUb(i) < minUcalc )
            {
                Xi[i] = this->getDomVarLb(i);
                Xj[i] = this->getDomVarLb(i);
            }

            else if (this->getDomVarLb(i) > maxUcalc)
            {
                Xi[i] = this->getDomVarUb(i);
                Xj[i] = this->getDomVarUb(i);
            }

            else // you assign them to get maximum deviation possible if nothing can be inferred
            {
               Xi[i] = this->getDomVarUb(i);
               Xj[i] = this->getDomVarLb(i);
            }
            sumXi += Xi[i] * Xi[i];
            sumXj += Xj[i];

        }
        sumXj = (sumXj * sumXj)/ this->getDomN();
        maxScalc = sumXi - sumXj; // n*variance
        maxScalc /= this->getDomN(); // variance
        maxScalc = sqrt(maxScalc); // std.deviation
        return maxScalc;
    }

    // This function creates the interval and it's values based on the current values of X
    void createInterval()
    {
        if(this->needDelete)
        {
//       cout << "Deleting previous intervals" << endl;
            for (int i = 0; i < this->numI; i++)
            {
                delete [] this->I[i]->M;
                delete [] this->I[i]->L;
                delete [] this->I[i]->R;
                delete this->I[i];
            }
            delete [] this->I;
        }
        int B[2*n]; // create an array to hold all current values of X
        for (int i = 0; i < n; i++)
        {
            B[i] = this->getDomVarLb(i);
            B[i+n] = this->getDomVarUb(i);
        }
        this->numI = this->mergesort(B, 2*n) - 1; // if B has 5 elements, this means you only have 4 intervals
        // Since now have all the unique domain variables, can create variable I
        this->I = new Interval*[this->numI];

        for (int i = 0; i < this->numI; i++)
        {
            // Set the upper and lower bound
            this->I[i] = new Interval(B[i], B[i+1]);
            // Now set values for M(Iq), R(Iq), L(Iq) as well as their cardinalities: m, r, l
            // While calculating ES
            this->I[i]->M = new int[this->getDomN()];
            this->I[i]->L = new int[this->getDomN()];
            this->I[i]->R = new int[this->getDomN()];
            this->I[i]->m = 0;
            this->I[i]->l = 0;
            this->I[i]->r = 0;
            this->I[i]->ES = 0;
            this->I[i]->C = 0;

            for (int j = 0; j <  this->getDomN(); j++)
            {
                if ( this->getDomVarUb(j) <= this->I[i]->Ilb )
                {
                    this->I[i]->L[this->I[i]->l] = j; // store the index of x which belongs to L(Iq)
                    this->I[i]->l++;
                    this->I[i]->ES += this->getDomVarUb(j);
                    this->I[i]->C += (this->getDomVarUb(j) * this->getDomVarUb(j));
                }
                else if (this->getDomVarLb(j) >= this->I[i]->Iub)
                {
                    this->I[i]->R[this->I[i]->r] = j; // store the index of x which belongs to R(Iq)
                    this->I[i]->r++;
                    this->I[i]->ES += this->getDomVarLb(j);
                    this->I[i]->C += (this->getDomVarLb(j) * this->getDomVarLb(j));
                }
                else
                {
                    this->I[i]->M[this->I[i]->m] = j; // store the index of x which belongs to M(Iq)
                    this->I[i]->m++;
                }
            }
            // Calculate the Vmin and Vmax for this interval
            // Vmin = ES + Imin*m
            this->I[i]->Vlb = this->I[i]->ES + (this->I[i]->Ilb * this->I[i]->m);
            this->I[i]->Vub = this->I[i]->ES + (this->I[i]->Iub * this->I[i]->m);

            // Calculate the range of n*variance that can be in this interval
            if (this->I[i]->m == 0)
            {
                // If m is 0, there means there are no Middle intervals
                this->I[i]->nVarLb = this->I[i]->C - ((this->I[i]->Vlb * this->I[i]->Vlb)/this->getDomN());
                this->I[i]->nVarUb = this->I[i]->C - ((this->I[i]->Vub * this->I[i]->Vub)/this->getDomN());
            }
            else
            {
                this->I[i]->nVarLb = this->I[i]->C + (((this->I[i]->Vlb - this->I[i]->ES)*(this->I[i]->Vlb - this->I[i]->ES))/this->I[i]->m) - ((this->I[i]->Vlb * this->I[i]->Vlb)/this->getDomN());
                this->I[i]->nVarUb = this->I[i]->C + (((this->I[i]->Vub - this->I[i]->ES)*(this->I[i]->Vub - this->I[i]->ES))/this->I[i]->m) - ((this->I[i]->Vub * this->I[i]->Vub)/this->getDomN());
            }

            // If special case where 2nd derivative is 0
            if ( this->I[i]->m == this->getDomN())
            {
                // do nothing, cause we're guaranteed the global min does not belong to this interval
                this->I[i]->q0 = -1; // m == n
            }
            else
            {
                // calculate the global min q if it was this interval
                this->I[i]->q0 = ((this->I[i]->ES * this->n)/(this->n  - this->I[i]->m)) ;
            }
            // Print only the first time
            if (this->needDelete == false)
            {
                // Print this to check interval initialize properly
                cout << "Interval: " << i << " Ilb: " << this->I[i]->Ilb << " Iub: " <<
                this->I[i]->Iub <<  " ES: " << this->I[i]->ES << " C: " << this->I[i]->C
                << " Vlb: " <<  this->I[i]->Vlb << " Vub: " << this->I[i]->Vub << endl
                << " nVarLb: " << this->I[i]->nVarLb << " nVarUb: " <<
                this->I[i]->nVarUb << " q0: " << this->I[i]->q0 << endl;
                cout << "m: " << this->I[i]->m << " l: " << this->I[i]->l << " r: " << this->I[i]->r <<endl;
            }
        }
        this->needDelete = true;
        return;
    }

    // Returns N*Variance value from the OPT graph given a value of q
    // can only be called after createInterval() function is used at least once
    // Note: It returns 0 for every q value outside the given domains of the intervals
    double OPT(double q)
    {
        double nVariance = 0;
        // If q is not within the intervals, return 0
        if((q < this->I[0]->Vlb) || (q > this->I[this->numI -1]->Vub))
        {
            return 0;
        }
        // Figure out which interval q is in
        int qIndex = 0;
        for (int i = 0; i < this->numI; i++)
        {
            if((q >= this->I[i]->Vlb) && (q <= this->I[i]->Vub))
            {
                qIndex = i;
                break;
            }
        }
        if (((double)this->I[qIndex]->m) == 0)
        {
            nVariance = (this->I[qIndex]->C - (pow((double) q, 2.0)/((double)this->getDomN()))) ;
        }
        else
        {
            nVariance = (this->I[qIndex]->C + ((pow((double)(q - this->I[qIndex]->ES), 2.0))/((double)this->I[qIndex]->m)) - (pow((double) q, 2.0)/((double)this->getDomN()))) ;
        }
        return nVariance;
    }

    // Note: optMax is from the SpreadConsS's upper bound
    // get the variance from S and multiply by N
    // bound = 0 => lower bound of interval
    // bound = 1 => upper bound of interval
    // Need bound cause for discontinuous graph,
    // bound can be upper or lower bound
    double DMAX(double q, double XMIN, double optMAX, int bound)
    {
        double dmax = 0;
        // If q is not within the intervals, return 0
        if((q < this->I[0]->Vlb) || (q > this->I[this->numI -1]->Vub))
        {
            return 0;
        }
        // Figure out which interval q is in
        int qIndex = 0;
        // If lower bound, start from lowest interval and look for it
        if (bound ==  0)
        {
            for (int i = 0; i < this->numI; i++)
            {
                if((q >= this->I[i]->Vlb) && (q <= this->I[i]->Vub))
                {
                    qIndex = i;
                    break;
                }
            }
        }
        // If upper bound, start from highest interval and look for it
        else // if (bound == 1)
        {
            for (int i = numI - 1; i >= 0; i--)
            {
                if((q >= this->I[i]->Vlb) && (q <= this->I[i]->Vub))
                {
                    qIndex = i;
                    break;
                }
            }
        }
        double v = 0;
        double m = 1.0*this->I[qIndex]->m;
        double a = 0;

        if (m == 0)
        {
            a = 1.0; // since the 1\m term will not exist if v is not involved
            v = 0; // v value is not needed since no X carries the value v
        }
        else
        {
            v = ((q - this->I[qIndex]->ES)/((double)m));
            a = 1.0 + (1.0/((double)m));
        }
        double n = this->getDomN();
        double opt = this->OPT(q); // takes into account if m is 0
        double b = XMIN - v;
        double c = opt - n*(pow((double) optMAX, 2.0));
        dmax = ((-b + (sqrt((b*b) - (a*c))))/((double)a));
        return dmax;
    }

    double DMIN(double q, double XMAX, double optMAX, int bound)
    {
        double dmin = 0;
        // If q is not within the intervals, return 0
        if((q < this->I[0]->Vlb) || (q > this->I[this->numI -1]->Vub))
        {
            return 0;
        }
        // Figure out which interval q is in
        int qIndex = 0;
        // If lower bound, start from lowest interval and look for it
        if (bound ==  0)
        {
            for (int i = 0; i < this->numI; i++)
            {
                if((q >= this->I[i]->Vlb) && (q <= this->I[i]->Vub))
                {
                    qIndex = i;
                    break;
                }
            }
        }
        // If upper bound, start from highest interval and look for it
        else // if (bound == 1)
        {
            for (int i = numI - 1; i >= 0; i--)
            {
                if((q >= this->I[i]->Vlb) && (q <= this->I[i]->Vub))
                {
                    qIndex = i;
                    break;
                }
            }
        }
        double m = 1.0*this->I[qIndex]->m;
        double v = 0;
        double a = 0;
        if (m == 0)
        {
            a  = 1.0;
            v  = 0;
        }
        else
        {
            v = ((q - this->I[qIndex]->ES)/((double)m));
            a = 1.0 + (1.0/((double)m));
        }
        double n = this->getDomN();
        double opt = this->OPT(q); // handles case where m == 0

        double b = v - XMAX; // note: This is different from DMAX
        double c = opt - n*(pow((double) optMAX, 2.0));
        dmin = ((-b + (sqrt((b*b) - (a*c))))/((double)a));
        return dmin;
    }

    double findDMax(double qSelected,int indexQSelected, int xIndex, int stdDevUb)
    {
   //     cout << " Parameters are: " << endl
   //     << "qSelected: " << qSelected << " indexQSelected: "
   //      << indexQSelected << " xIndex: " << xIndex << " stdDevUb: " << stdDevUb << endl;
    //    cout << " Calculated nVariance from stdDevUB: " << (pow((double) stdDevUb, 2.0)*this->n) <<endl;
        // Initialize d11 to maximum shift
        double d11 =  qSelected - this->I[indexQSelected]->Vlb;
        double dMax = 0;
        double a = 0;
        double b = 0;
        double c = 0;
        double v = 0;
        double nVarCurrent = 0;

        if (this->I[indexQSelected]->m == 0 )
        {
            v = 0;
            a = 1.0;
        }
        else
        {
            v = (qSelected - this->I[indexQSelected]->ES)/((double) this->I[indexQSelected]->m);
            nVarCurrent += (pow((double)(v - qSelected/((double)this->n)), 2.0) * this->I[indexQSelected]->m);
            a = 1.0+ 1.0/this->I[indexQSelected]->m;

        }

        int index = 0;
        for (int i = 0; i < this->I[indexQSelected]->l; i++)
        {
            index = this->I[indexQSelected]->L[i];
            nVarCurrent += pow((double)(this->getDomVarUb(index) - qSelected/((double)this->n)), 2.0);
        }
        for (int i = 0; i < this->I[indexQSelected]->r; i++)
        {
            index = this->I[indexQSelected]->R[i];
            nVarCurrent += pow((double)(this->getDomVarLb(index) - qSelected/((double)(this->n))),2.0) ;
        }

        b = this->getDomVarLb(xIndex) - v;
        c = nVarCurrent - (pow((double) stdDevUb, 2.0)*this->n);

    //    cout << " A: " << a << " B: " << b << " C: " << c <<" nVarCurrent " <<nVarCurrent << endl ;

        dMax = (-b + sqrt(pow (b, 2.0) - a*c))/a;
//cout << "dMax: " << dMax << " d11: " << d11 << endl;
        if (dMax < d11)
        {
            return dMax;
        }
        else
        {
            if((indexQSelected == 0) || (d11 == 0))
                return d11;
            else
            {
                // Update X's minimum value
                this->setDomVarLb(xIndex, getDomVarLb(xIndex) + floor(d11));
                this->createInterval(); // create the new intervals
                // Determine which indexQ is in
                for (int i = 0 ; i < this->numI ; i++)
                {
                    if ((qSelected >= this->I[i]->Vlb) && (qSelected <= this->I[i]->Vub))
                    {
                        indexQSelected = i;
                        break;
                    }
                }
                return (d11 + this->findDMax(qSelected, indexQSelected, xIndex, stdDevUb));
            }
        }
        // SHOULD NEVER COME HERE
        cout << "ERROR: WENT TO return dMAX" << endl;
        return dMax;
    }

    double findDMin(double qSelected,int indexQSelected, int xIndex, int stdDevUb)
    {
   //     cout << " Parameters are: " << endl
  //      << "qSelected: " << qSelected << " indexQSelected: "
    //     << indexQSelected << " xIndex: " << xIndex << " stdDevUb: " << stdDevUb << endl;
   //     cout << " Calculated nVariance from stdDevUB: " << (pow((double) stdDevUb, 2.0)*this->n) <<endl;

        // Initialize d11 to maximum shift
        double d11 =   this->I[indexQSelected]->Vub - qSelected;
        double dMin = 0;
        double a = 0;
        double b = 0;
        double c = 0;
        double v = 0;
        double nVarCurrent = 0;
        if(this->I[indexQSelected]->m == 0)
        {
            v = 0;
            a = 1.0;
        }
        else
        {
            v = (qSelected - this->I[indexQSelected]->ES)/((double) this->I[indexQSelected]->m);
            nVarCurrent += (pow((double)(v - qSelected/((double)this->n)), 2.0) * this->I[indexQSelected]->m);
            a = 1.0+ 1.0/this->I[indexQSelected]->m;

        }

        int index = 0;
        for (int i = 0; i < this->I[indexQSelected]->l; i++)
        {
            index = this->I[indexQSelected]->L[i];
            nVarCurrent += pow((double)(this->getDomVarUb(index) - qSelected/((double)this->n)), 2.0);
        }
        for (int i = 0; i < this->I[indexQSelected]->r; i++)
        {
            index = this->I[indexQSelected]->R[i];
            nVarCurrent += pow((double)(this->getDomVarLb(index) - qSelected/((double)(this->n))),2.0) ;
        }

        b = v - this->getDomVarUb(xIndex) ; // note: This is different from DMax
        c = nVarCurrent - (pow((double) stdDevUb, 2.0)*this->n);

    //    cout << " A: " << a << " B: " << b << " C: " << c <<" nVarCurrent " <<nVarCurrent << endl ;

        dMin = (-b + sqrt(pow (b, 2.0) - a*c))/a;
        cout << "dMin: " << dMin << " d11: " << d11 << endl;
        if (dMin < d11)
        {
            return dMin;
        }
        else
        {
            // return if last interval
            if((indexQSelected == (this->numI-1)) || (d11 == 0))
                return d11;
            else
            {
                // Update X's maximum value
                this->setDomVarUb(xIndex, getDomVarLb(xIndex) - floor(d11));
                this->createInterval(); // create the new intervals
                // Determine which indexQ is in
                for (int i = 0 ; i < this->numI ; i++)
                {
                    if ((qSelected >= this->I[i]->Vlb) && (qSelected <= this->I[i]->Vub))
                    {
                        indexQSelected = i;
                        break;
                    }
                }
                return (d11 + this->findDMin(qSelected, indexQSelected, xIndex, stdDevUb));
            }
        }
        // SHOULD NEVER COME HERE
        cout << "ERROR: WENT TO return dMin" << endl;
        return dMin;
    }

};

// End of Soons implementation set up without cons prop 
//-----------------------------------------------------------
// Thiago's implementation set up 
const double numerical_gap = 0.000001; // to avoid numerical methods error in SCIP according to Wen Yang





//#define pricing/delvars   	
/*
 * Data structure used in the spread constraint*/ 
 
// Create an interval data type 
struct interval
{
  double min, max, V_min, V_max, ES, GC;
  int *M;
  int *R;
  int *L;
  int m, r, l;       
};


// data for Spread constraints */
struct SCIP_ConsData
{
   // All your variables, its mean and the allowed maximum deviation 
   SCIP_VAR**     	 vars;         		//< array of variables
   int		 	 nvars;			//number of variables
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
// To handle bound change events 
   SCIP_EVENTHDLR*       eventhdlr;          //< event handler for bound change events 
// Do you use binary variables? 
   SCIP_Bool             usebinvars;         //< should the binary variables be used? 
   SCIP_Bool             edgefinding;        //< should edge-finding be executed? 
   SCIP_Bool             localcuts;          //< should cuts be added only locally? 
   SCIP_Bool             usecovercuts;       //< should covering cuts be added? 
   SCIP_Bool             sepaold;            //< shall old sepa algo be applied? 
   SCIP_Bool             fillbranchcands;    //< should branching candidates be added to storage? 
   SCIP_Bool             dualpresolve;       //< should dual presolving be applied? 
   SCIP_Bool             coeftightening;     //< should coeffisient tightening be applied? 
   SCIP_Bool             normalize;          //< should demands and capacity be normalized? */
   SCIP_Bool             usebdwidening;      //< should bound widening be used during conflict analysis? */
   SCIP_Longint          maxnodes;           //< number of branch-and-bound nodes to solve an independent disjunctive constraint  (-1: no limit) */
   SCIP_Bool             cutsasconss;        //< should the disjunctive constraint create cuts as knapsack constraints? */
   SCIP_Bool             coretimes;          //< should core-times be propagated (time tabling)? */
   SCIP_Bool             overload;           //< should edge finding be used to detect an overload? */
   SCIP_Bool             useadjustedjobs;    //< should during edge-finding jobs be adusted which run on the border of the effective time horizon? */
   SCIP_NODE*            lastenfolpnode;     /**< the node for which enforcement was called the last time (and some constraint was violated) */
   int                   nenfolprounds;      /**< counter on number of enforcement rounds for the current node */

};


//---------------------------------------------------------------------------------------------------------------------
/** Spread methods */
// This function calculates the maximum standard deviation that can result from a given set of domains of variables. 
// Used to reduce the upper bound of the spread constraint.
// returns the std. deviation calculated 
// This is used to prune the maximum std. deviation 
double max_stddev(int nvars, double vars[][2], double mean)
{ // Note: Nvars must be number of X variables 
   double sum_min = 0; // total sum of the lower bound
   double sum_max = 0; // total sum of the upper bound 
   for (int k=0;k<nvars;k++)
   {
      sum_min += vars[k][0]; // Smin
      sum_max += vars[k][1]; // Smax 
   }
   double min_mean = sum_min/nvars; // lower bound mean, Ulb
   double max_mean = sum_max/nvars; // upper bound mean, Uub 
   double x[nvars];		// Create an array of values for x
   double min_mean_;		// 
   double max_mean_;		// 
   double distance_x_max;	// 
   double distance_x_min;	// 

   // Assign each variable, x[k] to its extrema, such that the std dev is maximized
   // Assigning it to the extrema will maximize the std. dev as explained in the paper
   // "pg 13 of 15  Simplification and extension of Spread Constraint by Pierre et al" 
   // ERROR: THIS IS DEFINITELY WRONG CAUSE THE ALGORITHM TAKES O(n) instead of claimed by paper O(n^2) (refer to pg 13 of 15) 
   for (int k=0;k<nvars;k++)
   {
	  // assume x = xmax => min_mean increases 
      min_mean_ = min_mean + (vars[k][1] - vars[k][0])/nvars;

// As explained in the paper, find out which values to use before assigning values of x in order to get max std. deviation. 
		// lower bound of x = xmax
      distance_x_max =  std::min (fabs(vars[k][1] - min_mean_), // fabs(x) => absolute value of x 
	                         fabs(vars[k][1] - max_mean));  // min(x,y) => returns x if x <= y and y if y < x    

      //assume x = xmin => max_mean decreases 
      max_mean_ = max_mean - (vars[k][1] - vars[k][0])/nvars;
	  // upper bound of x = xmin 
      distance_x_min =  std::max (fabs(vars[k][0] - max_mean_), 
	                        fabs(vars[k][0] - min_mean));  
   
      if (distance_x_max > distance_x_min)
      {
			x[k] = vars[k][1]; // assign to upper bound extreme 
			min_mean = min_mean_; // update the minimum mean 
      }

	  // LINE HAHA SHOULD BE AN ERROR!!!!! cause suppose to check this instead
	  /*		
		// upper bound of x = xmax
		a = std::max (fabs(vars[k][1] - min_mean_), // fabs(x) => absolute value of x 
	                         fabs(vars[k][1] - max_mean));
		// lower bound of x = xmin
      b =  std::min (fabs(vars[k][0] - max_mean_), 
	                        fabs(vars[k][0] - min_mean));  		else if ( ) 
		else if (b > a) 
		{ 
		    x[k] = vars[k][0]; // assign to lower bound extreme 
			max_mean = max_mean_; // update the maximum mean
		} 
		
		else if (xmax < umin ) 
		=> assign x = xmin
		
		else if (xmin > umax)
		=> assign x = xmax 
		
		else 
		{
			// Nothing can be deduced from here. 
			calculate max_std_dev possible by assigning all xi and xj to whatever that was deduced
			and whatever that wasn't deduced, assign xi to its lower bound and xj to its upper bound 
		}
	  */ 
      else // LINE HAHA THIS SHOULDN'T BE CORRECT 
      {
         x[k] = vars[k][0]; // assign to lower bound extreme 
		 max_mean = max_mean_; // update the maximum mean
      }
	  // ERROR ENDS HERE! First, note that it does not check the opposite case where 
	  // 
   }
   double std_dev = 0;
   // Here, we're done assigning X to its appropriate extreme values, 
    // calculate the maximum std. deviatio that can occur 
   for(int i = 0; i < nvars; ++i )
      std_dev += pow(x[i] - mean, 2.0);  // numerator of variance calculation 
   std_dev = std_dev/nvars; // variance 
   std_dev = sqrt(std_dev); // std_deviation 
   return std_dev;
}
//---------------------------------------------------------------------------------------------------------------------
// This function is a simple bubble sort, but currently not correct as it does not take into account the last element. 
// helper function used by create_bounds()
// Yeap, the function it is used in may not really need to sort the last element. 
// Note: This bubblesort does not sort the last element cause it's only used in create_bounds method, 
// where the last element is guaranteed to be sorted. 
void sort(int n, double *a)
{
   bool sorted = false;
   double b;
   int current = n-1; 
   while (!sorted)
   {
      sorted = true;
      for (int i=0; i < current; i++)
         if (a[i] > a[i+1])
		{
	      b = a[i];	
		  a[i] = a[i+1];
		  a[i+1] = b;
		  sorted = false;
        }
	  current--;
   }
}
//---------------------------------------------------------------------------------------------------------------------
// Given the domains of the variables, creates an ordered set of the existing values
// vars => variables with their lower bounds and upper bounds
// nb => Number of bounds
// Basically, from all x1, x2, x3, it returns the new I1, I2, I3, which can be all the way up to I(2n-1) 
// Creates the Contiguous Interval as explained in the (paper pg 3 of 15)
// Note: This functio works correctly as checked (Soon) 

void create_bounds(int n, double vars[][2], int &nb, double* &bounds)
{
   int l, k;
   bool new_min, new_max;
   double min_value = vars[0][0];
   double max_value = vars[0][1];
   nb = 0;   
// Get min_value and max_value 
   for (l = 0; l < n; l++) 
   {     
      // Update min and max value 
      if (vars[l][0] < min_value)  min_value = vars[l][0];
      if (vars[l][1] > max_value)  max_value = vars[l][1];

      new_min = true;
      new_max = true; 
      for (k = 0; k < nb; k++)
      {
         if (vars[l][0] == bounds[k]) 
	 {
	    new_min = false;
         }
         if (vars[l][1] == bounds[k])
	 {
 	    new_max = false;
         }
      }
      if (new_min)
      {
         bounds[nb] = vars[l][0];  
         nb++;
      }
      if (new_max)
      {
         bounds[nb] = vars[l][1];  
         nb++;
      }
   }
printf("PRINTING BOUNDS before Sort\n"); 
for (k = 0; k< nb; k++)
{
printf("Bound %d is %f \n", k, bounds[k]); 
}


   sort(nb, bounds);
printf("PRINTING BOUNDS after Sort\n"); 
for (k = 0; k< nb; k++)
{
printf("Bound %d is %f \n", k, bounds[k]); 
}
}

//---------------------------------------------------------------------------------------------------------------------
// This function returns indices V, ES, R, L, M... for a given interval I, Note: 
// If you have 4 different interval I as shown in (pg 5 of 15), you will need to call this function 4 times, 
// once for each interval
// Function used in the propagator
// Takes in number of variables, each with lower and upper bound and a pointer for the intervals 

void calculate_indices(int nvars, double vars[][2], interval &I)
{ 
   //A = the values assumed by variables on each interval I (not in order)
   I.l=0; I.r=0; I.m=0; // Actual number of variables for left, right, and middle. 
   I.L = new int[nvars]; // Initialize space for nvars although I.L can be <= nvars 
   I.M = new int[nvars];
   I.R = new int[nvars];     
   I.ES = 0;            
   for (int i = 0; i < nvars; i++) 
   {
      // If upper bound is lower than minimum, it must be the L(I) 
      if (vars[i][1] <= I.min)
      {
	     I.L[I.l] = i;
	     I.l++;
	     I.ES += vars[i][1];
      }
      // Else if lower bound is higher than maximum, it must be R(I) 
      else if (vars[i][0] >= I.max)
      {	 
	     I.R[I.r] = i;
	     I.r++;
	     I.ES += vars[i][0];
      }
      else // else, it must be M(I) 
      {
  	     I.M[I.m] = i;
         I.m++; // cardinality of M increases
      }
   }
   I.V_min = I.ES + I.min*I.m;  // the Vlowerbound for this interval
   I.V_max = I.ES + I.max*I.m;  // the Vupperbound for this interval 
   I.GC = (double) nvars*I.ES/((double)nvars - I.m);
}

//---------------------------------------------------------------------------------------------------------------------
// This function calculates the optimal value for std. deviation for the graph after accepting a specific value of q. 
//calculates optPI (std deviation) for a given interval and a given q value (q = nvars * mean).  Function used in the propagator
// Basically calculates n*variance = (xi-u)^2 where xi is either max, min or v. 
double Calculate_optPI(interval Iq, double q, double vars[][2], int nvars, double &v)
{
   double optPI = 0;
   for (int i = 0; i < Iq.l; i++) 
     optPI += pow( (double)vars[Iq.L[i]][1] - (double) q/(double)nvars, 2.0);
   //sums variables that belong to R(I)
   for (int i = 0; i < Iq.r; i++) 
     optPI += pow( (double)vars[Iq.R[i]][0] - (double) q/(double)nvars, 2.0);
   //sums variables that belong to M(I).      m.v =  q - ES
   if (Iq.m > 0)  
   {
     v = (q - Iq.ES)/((double)Iq.m);
     optPI += pow( v - q/(double)nvars, 2.0) * Iq.m;
   }
   return optPI;
}

//-------------------------------------------  --------------------------------------------------------------------------
// NOWNOWNOWNOWNOW : Understand exactly what is being done and analytically analyze if this makes sense  
//Finds the maximum shift (d) that can be applied to a variable's domain such that PiMax (maximum std dev allowed) is reached. Used to update upper bound of the variable.
// The algorithm is illustrated and explained in (pg 6 of 15) 
// Although you do not fully understand the recursive call backwards yet. 
double FindDMax(int nvars, double vars[][2], int index_var, int index_int, double pi1Max, double q, interval Iq){
//   printf("\nFinddMax");
   double dMax;
   double n = (double) nvars;
   double v;
   double optPI = Calculate_optPI(Iq, q, vars, nvars, v);
   double a, b, c, d1;
   d1 = q - Iq.V_min;
   a = 1.0 + 1.0/Iq.m;
   b = vars[index_var][0] - v;
   c = optPI - pi1Max;

   dMax = (-b + sqrt(pow (b, 2.0) - a*c))/a;
   if (dMax < d1)
	return dMax;
   else
   {
    	if (index_int == 0) 
    	    return d1;
    	else
    	{

	    if (d1 <= -numerical_gap)
		return INT_MAX;

    	    //calculate new v and new optPI
    	    vars[index_var][1] = vars[index_var][1] + d1;
    	    vars[index_var][0] = vars[index_var][0] + d1;
    	    
    	    //creates new problem after adding new variable
            int nb_;
            double *bounds_ = new double[nvars*2+5]; 
            create_bounds(nvars, vars, nb_, bounds_);
            interval* I_ = new interval[nb_-1+5];
            int j_;
	    int nI_ = nb_-1;
	    bool int_found = false;	
            for (int k=0;k < nb_-1; k++) //creates new intervals, since the domains of a variable changed
            {
                I_[k].min = bounds_[k];
                I_[k].max = bounds_[k+1];	
                if ((bounds_[k] == Iq.min) && (bounds_[k+1] == Iq.max) && (!int_found))	 //the same Iq was found, so we store its index
		{
		   int_found = true; 
                   j_ = k;
		}
            }

//the interval changed when the variable was shifted from x to x'.. so the old Iq doesn't exist anymore. We are now looking for the Iq with its updated lower bound.
	    if (!int_found) 
	    {
                for (int k=0;k < nb_-1; k++)
                    if ((bounds_[k] == Iq.min) && (bounds_[k+1] <= vars[index_var][0]))
                    {
			j_ = k;
			break;
 		    }
	    }

	    if (j_ == 0) //first interval.. will result in Segmentation fault since there's no predecessor
   	       index_int = j_+1;  	
  	    else
  	       index_int = j_; 

	    //prevent segment fault errors
	    if ((index_int-1 >= nI_) || (index_int-1 <=-1))
		return INT_MAX;

            calculate_indices(nvars, vars, I_[index_int-1]);
  // NOWNOWNOWNOWNOW : Understand what past student implemented and check if it makes sense. 
	    d1 = d1 + FindDMax(nvars, vars, index_var, index_int-1, pi1Max, q, I_[index_int-1]);
	    delete [] bounds_;
  	    delete [] I_[index_int-1].L;
	    delete [] I_[index_int-1].R;
	    delete [] I_[index_int-1].M;
	    delete [] I_;

    	    return d1;
    	}
   }
}
//---------------------------------------------------------------------------------------------------------------------
//Finds the minimum shift (-d) that can be applied to a variable's domain such that PiMax (maximum std dev allowed) is reached. Used to update lower bound of the variable.
// This algorithm is similar to FindDMax as the problem is similar, so just need to reimplemented in the other direction. 
// Note: Have to manually derive this algorithm from findDMax by changing some equality operators and variables. 
// Similar to (pg 6 of 15)
double FindDMin(int nvars, double vars[][2], int index_var, int index_int, double pi1Max, double q, interval Iq, int last_interval){
   double dMin;
   double n = (double) nvars;
   double v;
   double optPI = Calculate_optPI(Iq, q, vars, nvars, v);
   double a, b, c, d1;
   d1 = - q + Iq.V_max;
   if (Iq.m > 0)
   {
      a = 1.0 + 1.0/Iq.m;
      b = - vars[index_var][1] + v;
      c = optPI - pi1Max;
      dMin = (-b + sqrt(pow (b, 2.0) - a*c))/a;
   }
   else
      dMin = INT_MAX;   

   if (dMin < d1)
	return dMin;
   else
   {
    	if (index_int == last_interval) //last interval??
    	    return d1;
    	else
    	{

	    if (d1 <= -numerical_gap)
		return INT_MAX;

    	    //calculate new v and new optPI
    	    vars[index_var][1] = vars[index_var][1] - d1;
    	    vars[index_var][0] = vars[index_var][0] - d1;
    	    
    	    //creates new problem after adding new variable
            int nb_;
            double *bounds_ = new double[nvars*2+5];
            create_bounds(nvars, vars, nb_, bounds_);
            interval* I_ = new interval[nb_-1+5];
            int j_;
	    bool int_found = false;		
            int nI_ = nb_ - 1;

            for (int k=0;k < nI_; k++) //creates new intervals, since the domains of a variable changed
            {
                I_[k].min = bounds_[k];
                I_[k].max = bounds_[k+1];
                if ((bounds_[k+1] == Iq.max) && (bounds_[k] == Iq.min) && (!int_found))  //the same Iq was found, so we store its index
		{
		   int_found = true;
                   j_ = k;
		}
            }


//the interval changed when the variable was shifted from x to x'.. so the old Iq doesn't exist anymore. We are now looking for the Iq with its updated lower bound.
	    if (!int_found) 
	    {
                for (int k=0;k < nI_; k++)
                    if ((bounds_[k+1] == Iq.max) && (bounds_[k] >= vars[index_var][1]))
                    {
			j_ = k;
			break;
 		    }
	    }

	    if (j_ == nI_-1) //last interval.. will result in Segmentation fault since there's no successor
   	       index_int = j_-1;  	
  	    else
  	       index_int = j_; 


	    //prevent segment fault errors
	    if ((index_int+1 >= nI_) || (index_int+1 <=-1))
		return INT_MAX;

            calculate_indices(nvars, vars, I_[index_int + 1]);

            d1 = d1 + FindDMin(nvars, vars, index_var, index_int+1, pi1Max, q, I_[index_int+1], nI_-1);
	    delete [] bounds_;
  	    delete [] I_[index_int+1].L;
	    delete [] I_[index_int+1].R;
	    delete [] I_[index_int+1].M;
	    delete [] I_;
	    return d1;
    	}
   }
}

//---------------------------------------------------------------------------------------------------------------------
/** Local methods */

static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to consdata */
   SCIP_VAR**            vars,               /**< array of integer variables */
   int 			 nvars
   )
{

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(vars != NULL);	
   assert(nvars > 0);
   assert(mean != NULL);
   assert(max_deviation != NULL);

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   
   /* Marks that all variables must not be multi-aggregated - since Spread is a global constraint that tightens more effectively when more variables are included in the constraint. */ 
   if (nvars > 0){
   	for (int i=0; i < nvars; i++){
	   SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, vars[i]) );
        }
   }
//   SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, max_deviation) );

   (*consdata)->nvars = nvars;
   (*consdata)->vars = vars;
   
   if( nvars > 0 )
   {
      assert(vars != NULL); /* for flexelint */

      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );

      /* transform variables, if they are not yet transformed */
      if( SCIPisTransformed(scip) )
      {
         SCIPdebugMessage("get tranformed variables and constraints\n");

         /* get transformed variables and do NOT capture these */
         SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );
      }
   }
   else
   {
printf("Consdata is null\n"); 
      (*consdata)->vars = NULL;
   }

   return SCIP_OKAY;
}
//---------------------------------------------------------------------------------------------------------------------

// QUESTION: NO IDEA IF NEED THIS 
static
SCIP_RETCODE Relaxation_Mean(SCIP* scip, int nvars, SCIP_VAR** vars, double mean, const char* name)
{
   // Adds relaxations for the spread constraint
   SCIP_CONS* cons_mean;
   char cons_name[100];
   snprintf(cons_name, sizeof name, "Relaxation_%s", name);  
   SCIP_Real value = mean*nvars;
   SCIP_CALL( SCIPcreateConsLinear(scip, & cons_mean, name, 0, NULL, NULL, value, value, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   for (int j=0; j<nvars;j++)
       SCIP_CALL( SCIPaddCoefLinear(scip, cons_mean, vars[j], 1.0) );  
   SCIP_CALL( SCIPaddCons(scip, cons_mean) );
   return SCIP_OKAY;
}

//---------------------------------------------------------------------------------------------------------------------

/** creates constaint handler data for Spread constraint handler */
static
SCIP_RETCODE conshdlrdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata,       /**< pointer to store the constraint handler data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   /* create precedence constraint handler data */
   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPallocMemory(scip, conshdlrdata) );

   /* set event handler for checking if bounds of start time variables are tighten */
   (*conshdlrdata)->eventhdlr = eventhdlr;

#ifdef SCIP_STATISTIC
   (*conshdlrdata)->nirrelevantjobs = 0;
   (*conshdlrdata)->nalwaysruns = 0;
   (*conshdlrdata)->ndualfixs = 0;
   (*conshdlrdata)->nremovedlocks = 0;
   (*conshdlrdata)->ndecomps = 0;
   (*conshdlrdata)->nallconsdualfixs = 0;
#endif

   return SCIP_OKAY;
}


//---------------------------------------------------------------------------------------------------------------------
static
SCIP_RETCODE consdataCatchEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< disjunctive constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   int v;
   /* catch event for every single variable */
   for( v = 0; v < consdata->nvars; ++v )
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[v],
            SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );

							// std Dev. 
 /*  SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[consdata->nvars-2],
            SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );*/
   return SCIP_OKAY;
}

//---------------------------------------------------------------------------------------------------------------------
SCIP_RETCODE SCIPcreateConsSpread(
   SCIP* 		 scip, /**< SCIP data structure */
   SCIP_VAR** 		 vars,                /**< variables */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
)
{ 
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   assert(scip != NULL);
   assert(vars != NULL);
   assert(cons != NULL);
   assert(nvars > 0);
   assert(mean != NULL);
   assert(max_deviation != NULL);

   /* find the Spread constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage(""CONSHDLR_NAME" constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   SCIPdebugMessage("create Spread constraint <%s> with %d variables\n", name, nvars);


   int i;	
   for(i=0;i<nvars;i++)
   {
      printf("Spread constraint var %s with [%f, %f] \n", SCIPvarGetName(vars[i]), SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i]));
   }

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, vars, nvars));

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata,
         initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   /* create relaxation for the mean */
//   Relaxation_Mean(scip, nvars, vars, mean, name); QUESTION: DO I NEED THIS? 
   if( SCIPgetStage(scip) != SCIP_STAGE_PROBLEM )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      /* get event handler */
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* catch bound change events of variables */
      SCIP_CALL( consdataCatchEvents(scip, consdata, conshdlrdata->eventhdlr) );
   }

   return SCIP_OKAY;
}

//---------------------------------------------------------------------------------------------------------------------
/*
 * Callback methods of constraint handler
 */

/* TODO: Implement all necessary constraint handler methods. The methods with #if 0 ... #else #define ... are optional */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopySpread)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrSpread(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}
#else
#define conshdlrCopySpread NULL
#endif
//---------------------------------------------------------------------------------------------------------------------

/** frees constraint handler data for logic or constraint handler */
static
void conshdlrdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

   SCIPfreeMemory(scip, conshdlrdata);
}

//---------------------------------------------------------------------------------------------------------------------
/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeSpread)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdataFree(scip, &conshdlrdata);

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
#if 0
static
SCIP_DECL_CONSINIT(consInitSpread)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Spread constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitSpread NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_CONSEXIT(consExitSpread)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Spread constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitSpread NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreSpread)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Spread constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
   printf("\n\nCONSInitPre");
   return SCIP_OKAY;
}
#else
#define consInitpreSpread NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitpreSpread)
{  
   printf("\n\nCONSExitPre");

   return SCIP_OKAY;
}
#else
#define consExitpreSpread NULL
#endif

//---------------------------------------------------------------------------------------------------------------------
/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
//#if 0
static
SCIP_DECL_CONSINITSOL(consInitsolSpread)
{  /*lint --e{715}*/
/*   SCIPerrorMessage("method of Spread constraint handler not implemented yet\n");
   SCIPABORT(); */ /*lint --e{527}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->lastenfolpnode = NULL;
   conshdlrdata->nenfolprounds = 0;
   return SCIP_OKAY;
}
/*#else
#define consInitsolSpread NULL
#endif*/


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_CONSEXITSOL(consExitsolSpread)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* release the rows of all constraints */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* free rows */
   //   SCIP_CALL( consdataFreeRows(scip, &consdata) );
   }

   return SCIP_OKAY;
}
#else
#define consExitsolSpread NULL
#endif


/** frees specific constraint data */
#if 0
static
SCIP_DECL_CONSDELETE(consDeleteSpread)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Spread constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeleteSpread NULL
#endif


/** transforms constraint data into data belonging to the transformed problem */
#if 0
static
SCIP_DECL_CONSTRANS(consTransSpread)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Spread constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransSpread NULL
#endif


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
#if 0
static
SCIP_DECL_CONSINITLP(consInitlpSpread)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Spread constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpSpread NULL
#endif

//---------------------------------------------------------------------------------------------------------------------
// separation method of constraint handler for LP solutions 
static
SCIP_DECL_CONSSEPALP(consSepalpSpread)
{  

  return SCIP_OKAY;

}


/** separation method of constraint handler for arbitrary primal solutions */
#if 0
static
SCIP_DECL_CONSSEPASOL(consSepasolSpread)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Spread constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepasolSpread NULL
#endif
//---------------------------------------------------------------------------------------------------------------------
//function called by CheckConsSpread
//Checks if sum of the values of the variables divided by n results in the proposed mean (failure type 1)
//Checks if std dev is within proposed gap (failure type 2)

static SCIP_RETCODE CheckCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to be checked, or NULL for current solution */
   int                   &violated	            /**< pointer to store whether the constraint is violated */
   )
{ 
   assert( consdata != NULL );
   violated = 0;	
   SCIP_VAR** vars;
   vars = consdata->vars;
   assert( consdata->nvars != NULL );	
   double mean = 0;

	    printf("Bounds at CheckCons are\n");
            for (int k = 0; k < consdata->nvars-2; k++)
            {
		printf("X[ %d ] LB= %f  UB= %f \n", k, SCIPvarGetLbLocal(consdata->vars[k]), SCIPvarGetUbLocal(consdata->vars[k])); 
	    }
	    printf("STDDEV LB= %f  UB= %f \n", SCIPvarGetLbLocal(consdata->vars[consdata->nvars-2]), SCIPvarGetUbLocal(consdata->vars[consdata->nvars-2])); 
	    printf("MEAN LB= %f  UB= %f \n", SCIPvarGetLbLocal(consdata->vars[consdata->nvars-1]), SCIPvarGetUbLocal(consdata->vars[consdata->nvars-1])); 

   for(int i = 0; i < ((int)consdata->nvars-2); i++ ) // First nvars-2 values are X, nvars-1 is stdDev, nvars-1 is mean 
   {
   	mean += SCIPgetSolVal(scip, sol, vars[i]);
	printf("Solution value for each X[ %d ] is %f\n", i, SCIPgetSolVal(scip, sol, vars[i])); 
   }

					// int W = consdata->W
					// SCIPgetSolVal(scip, sol, W) ; 
// Calculate the final mean from the solution values of all the x 
   mean = mean/((double)consdata->nvars-2);


// If the mean is not the same as the values of mean assigned, it is violated, 
// return error. 
//mean consistency

printf("Solution value for Mean is %f \n", SCIPgetSolVal(scip,sol,consdata->vars[consdata->nvars-1])); 
   if (mean != SCIPgetSolVal(scip,sol,consdata->vars[consdata->nvars-1]))
   {
	violated = 1;
        return SCIP_OKAY;
   }


// Now calculate the std. deviation 
//max dev consistency
   double std_dev = 0;
   for(int i = 0; i < ((int)consdata->nvars-2); ++i )
	std_dev += pow(SCIPgetSolVal(scip, sol, vars[i]) - mean, 2.0);
   std_dev = std_dev/(consdata->nvars-2);
   std_dev = sqrt(std_dev);
   double upper_bound = SCIPvarGetUbLocal(consdata->vars[consdata->nvars-2]); // get upper bound of std. deviation 

// If the std. deviation is larger than the upper bound of the std. deviation, it is violated. 
   if (std_dev > upper_bound + numerical_gap*10)
	violated = 2;

// Note: Does not account for any min_std deviation, max_deviation is good enough. 
   return SCIP_OKAY;
}

// constraint enforcing method of constraint handler for LP solutions 
#if 0
static
SCIP_DECL_CONSENFOLP(consEnfolpSpread)
{  


   return SCIP_OKAY;

}
#else
#define consEnfolpSpread NULL
#endif

// constraint enforcing method of constraint handler for pseudo solutions 
#if 0
static
SCIP_DECL_CONSENFOPS(consEnfopsSpread)
{  


   return SCIP_OKAY;
}
#else
#define consEnfopsSpread NULL
#endif

// feasibility check method of constraint handler for solutions 
static
SCIP_DECL_CONSCHECK(consCheckSpread)
{  
   int violated = 0;
   int c, i, j;
   SCIPdebugMessage("\nConsCheck");   

   assert( conss != NULL );
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   // Get all the Spread Constraint that exist 
   for (c = 0; c < nconss; ++c)
   {
	 SCIP_CONSDATA* consdata;
         consdata = SCIPconsGetData(conss[c]);
         SCIP_CALL ( CheckCons(scip, consdata, sol, violated) );
	 if (violated != 0)  break;
   }
   if( violated != 0 )
   {
      *result = SCIP_INFEASIBLE;
      if (violated == 1) //violated mean
        SCIPinfoMessage(scip, NULL, "Violated mean of Spread constraint");
      else
         SCIPinfoMessage(scip, NULL, "Violated standard deviation of Spread constraint");
      
   }
   else
      *result = SCIP_FEASIBLE;
   return SCIP_OKAY;

}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropSpread)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);
   SCIP_Bool infeasible;
   (*result) = SCIP_DIDNOTFIND;

printf("ENTERING CONSPROP\n"); 

// loop through all Spread Constraints  
   for (int c = 0; c < nconss; ++c)
   {

   double minSTDDEV = 99999;
 double tempMean = 0;
 double haha = 0;
double tempX1 = 0;
double tempX2 = 0;
double tempX3 = 0;
    // Temporary to brute force to calculate optimal answer
    for(int i = 1; i < 4; i++)
    {
        for(int j = 2; j<7; j++)
        {
            for(int k = 3; k<9; k++)
            {
                tempMean = (i+j+k)/3.0;
                if(tempMean == 5)
                {
                    haha = ((i - tempMean)*(i - tempMean) +  (j - tempMean)*(j - tempMean) +  (k - tempMean)*(k - tempMean));
                    haha /= 3.0;
                    haha = sqrt(haha);
                    if (haha < minSTDDEV)
                    {
                        minSTDDEV = haha;
                        tempX1 = i;
                        tempX2 = j;
                        tempX3 = k;
                    }
                }
            }
        }
    }
    cout << " Optimal value for X1: " << tempX1 << " X2: " << tempX2 << " X3: " << tempX3 << endl;

    // TODO Reverse Engineer: Brute Force find optimal solution for X1, X2, X3 without U and S
    // Find the minimum S able to be obtained,
    // and check if propagation prunes the optimal solution.

    cout << "Hello world!" << endl;
    int i = 0;

// Currently, 12th August 2014,
// my propagation algorithm appears to be correct (from brute force calculations)
// and prunes tighter than CP Optimizer for both general and specific mean case


// Assuming read N, lb and ub values from SCIP
/*
    int N = 3; int lb[3] ={1, 2, 3}; int ub[3] = {3,6,9};
    domainVar SpreadConsX(N, 'X'); // create space for 3 domain variables
    for (int i = 0; i < N; i++)
    {
        SpreadConsX.updateDomainVar(i, lb[i], ub[i]);
    }
	*/
	        SCIP_CONSDATA* consdata;
		SCIP_CONS* cons;
     	cons = conss[c];
     	assert( cons != NULL );	
     	consdata = SCIPconsGetData(cons);
		// Copy paste from SCIP into propagation 
	// Get number of variables 
        int N = (consdata->nvars-2); // number of Variables for X 
		domainVar SpreadConsX(N, 'X'); // create space for 3 domain variables
		for (int i = 0; i < N; i++)
		{
			SpreadConsX.updateDomainVar(i, SCIPvarGetLbLocal(consdata->vars[i]), SCIPvarGetUbLocal(consdata->vars[i]));
		}		
		SpreadConsX.printDomVar();
		domainVar SpreadConsU('U', 1, SCIPvarGetLbLocal(consdata->vars[consdata->nvars-2]), SCIPvarGetUbLocal(consdata->vars[consdata->nvars-2]));
		SpreadConsU.printDomVar();
		domainVar SpreadConsS('S', 1, SCIPvarGetLbLocal(consdata->vars[consdata->nvars-1]), SCIPvarGetUbLocal(consdata->vars[consdata->nvars-1]));
		SpreadConsS.printDomVar();
		// Create the intervals
		SpreadConsX.createInterval();
		// Done getting domain variables from SCIP, now start solving 
		// Note: I did not handle cases where no feasible solution is found to output to SCIP 

/*
    domainVar SpreadConsX('X', 3, 10, 20, 12, 20, 3, 30);
    SpreadConsX.printDomVar();
    domainVar SpreadConsU('U', 1, 1, 9);
    SpreadConsU.printDomVar();
    domainVar SpreadConsS('S', 1, 1, 30);
    SpreadConsS.printDomVar();
    // Create the intervals
    SpreadConsX.createInterval();
*/

/*
    domainVar SpreadConsX('X', 3, 1, 3, 2, 6, 3, 9);
    SpreadConsX.printDomVar();
    domainVar SpreadConsU('U', 1, 2, 2);
    SpreadConsU.printDomVar();
    domainVar SpreadConsS('S', 1, 1, 10);
    SpreadConsS.printDomVar();
    // Create the intervals
    SpreadConsX.createInterval();
*/
/*
// Tested and appears to be correct
// For debugging OPT
    for (int yp = 0; yp < 19; yp++)
    {
        double hehe = SpreadConsX.OPT(yp);
        cout << "q: " << yp << " OPT: " << hehe << endl;
    }
*/
/*
// Tested and appears to be correct
// For debugging DMAX Graph
    for (int hp = 0; hp < 21; hp++)
    {
        // check lower bound
        double huhuLB = SpreadConsX.DMAX(hp, 1.0, 50, 0);
        cout << "q: " << hp << " OPTlb: " << huhuLB << endl;
        double huhuUB = SpreadConsX.DMAX(hp, 1.0, 50, 1);
        cout << "q: " << hp << " OPTub: " << huhuUB << endl;
    }
*/

/*
    domainVar SpreadConsX('X', 3, 1, 3, 2, 6, 3, 9);
    SpreadConsX.printDomVar();
    domainVar SpreadConsU('U', 1, 5,5);
    SpreadConsU.printDomVar();
    domainVar SpreadConsS('S', 1, 1, 100);
    SpreadConsS.printDomVar();
    // Create the intervals
    SpreadConsX.createInterval();
*/

    // Note: Found out need to prune both ways for M(Iq)
    // Note: Found out can only change bounds after pruning all X, so that don't result in
    // QIndex not existing segmentation fault.
    // Also found out that adding the second way to prune M(Iq) can result in
    // S being (0,0) as it prunes further and when that happens,
    // I cannot run any X propagation anymore, so break right away
    // Note: Found out that forgot to handle special case where m = 0,
    // if m = 0, cannot divide by m!! thus, that's why previous code got lazy and just added 0.000001
    // but the right thing is to handle everything where m = 0 => need to recalculate a simpler equation
    // which was done and now it works.

/*
    domainVar SpreadConsX('X', 3, 1, 3, 2, 6, 3, 9);
    SpreadConsX.printDomVar();
    domainVar SpreadConsU('U', 1, 6, 6);
    SpreadConsU.printDomVar();
    domainVar SpreadConsS('S', 1, 1, 4);
    SpreadConsS.printDomVar();
    // Creation of intervals is where R, M, L is determined for each of the original intervals
    // Note: Creation of intervals only depends on X and not S or U
    // Create the intervals
    SpreadConsX.createInterval();
*/

    // Propagation Algorithm
    bool changed = true;
    // While something is pruned , keep running propagation
    while (changed)
    {
        // Initialize to false
        changed = false;

        //------------------------------------------------------------
        // Prune the Std Deviation Maximum value using X
        double maxScalc  = 0;
        if ( SpreadConsU.getDomVarLb(0) == SpreadConsU.getDomVarUb(0))
        {
            // Specifc case, faster but not as tight
            maxScalc = SpreadConsX.maxStdDevSpecific( (double) SpreadConsU.getDomVarUb(0));
            // General case, prunes tighter
            maxScalc = SpreadConsX.maxStdDev();
        }
        // General Case
        else
        {
            maxScalc = SpreadConsX.maxStdDev();
           //  maxScalc = SpreadConsX.maxStdDevThiago();// Note: Has been confirmed to be wrong
        }
        if (maxScalc < SpreadConsS.getDomVarUb(0))
        {
           // changed = true; // Note: Dont need to re-run since this is the first to get propagated
            SpreadConsS.setDomVarUb(0,floor(maxScalc));
            if (SpreadConsS.getDomVarUb(0) < SpreadConsS.getDomVarLb(0))
            {
                SpreadConsS.setDomVarLb(0,SpreadConsS.getDomVarUb(0));
            }
        }
        cout << "Propagate S using X" <<endl;
        SpreadConsS.printDomVar();
        //------------------------------------------------------------
        // Prune the Mean bounds using X
        double Sub = 0; // Sum of upper bounds
        double Slb = 0; // Sum of lower bounds

        for (i = 0; i < SpreadConsX.getDomN(); i++ )
        {
            Slb += SpreadConsX.getDomVarLb(i);
            Sub += SpreadConsX.getDomVarUb(i);
        }
        double minUcalc = Slb/SpreadConsX.getDomN();
        double maxUcalc = Sub/SpreadConsX.getDomN();
        if (minUcalc > SpreadConsU.getDomVarLb(0))
        {
            if (minUcalc > SpreadConsU.getDomVarUb(0))
            {
                cout << "No feasible solution as minUCalc from X is > Umax"  << endl;
                return 0;
            }
            SpreadConsU.setDomVarLb(0, ceil(minUcalc));
            changed = true;
        }
        if (maxUcalc < SpreadConsU.getDomVarUb(0))
        {
            if (maxUcalc < SpreadConsU.getDomVarLb(0))
            {
                cout << "No feasible solution as maxUCalc from X is < Umin"  << endl;
                return 0;
            }
            SpreadConsU.setDomVarUb(0, floor(maxUcalc));
            changed = true;
        }
        cout << "Propagate U using X" <<endl;
        SpreadConsU.printDomVar();
        //------------------------------------------------------------
        // Prune the Mean Bounds using S

        // Calculate the max N*Variance based on given std. deviation S
        // Note: For SCIP , due to Bound Consistency, can only prune using max variance and not minimum.
        double nVarUbCalc = SpreadConsX.getDomN() * SpreadConsS.getDomVarUb(0) * SpreadConsS.getDomVarUb(0);
//   cout << "nVarUbCalc " << nVarUbCalc << endl;

        // Check for consistency of solution
        double lowestPossibleNVar = SpreadConsX.I[0]->nVarLb; // Initialize to first lower bound
        for (int i = 0; i < SpreadConsX.numI; i++)
        {
            // If this interval contains the global minimum
            if ((SpreadConsX.I[i]->q0 <= SpreadConsX.I[i]->Vub) && (SpreadConsX.I[i]->q0 >= SpreadConsX.I[i]->Vlb))
            {
                lowestPossibleNVar = SpreadConsX.I[i]->C + (((SpreadConsX.I[i]->q0 - SpreadConsX.I[i]->ES)*(SpreadConsX.I[i]->q0 - SpreadConsX.I[i]->ES))/SpreadConsX.I[i]->m) - ((SpreadConsX.I[i]->q0 * SpreadConsX.I[i]->q0)/SpreadConsX.getDomN());
                break; // break out of for loop cause guaranteed to be minimal
            }
            // If it does not contain the global minimum, one of the end points must be reaching the global minimum,
            // thus, set it to one of the end points.
            else
            {
                if(SpreadConsX.I[i]->nVarLb < lowestPossibleNVar)
                {
                    lowestPossibleNVar = SpreadConsX.I[i]->nVarLb;
                }
                if(SpreadConsX.I[i]->nVarUb < lowestPossibleNVar)
                {
                    lowestPossibleNVar = SpreadConsX.I[i]->nVarUb;
                }
            }
        }
        cout << "nVarUBCalc is " << nVarUbCalc << " whereas lowestPossibleNvar is " << lowestPossibleNVar << endl;
        if (lowestPossibleNVar > nVarUbCalc)
        {
        // Note: Can't do this as it turns out the algorithm below doesn't work
        // on the specific example of [1,3], [2,6]  , [3,9] variables
        // It will end up pruning to [3,3] for all 3 variables
        // which prunes S to [0,0], which makes algorithm doesn't work.
        // However, still need above code to solve for nVarUbCalc
        // which is used below
        //    cout << "NO FEASIBLE SOLUTION " << endl
        //       << "as maximum n*Variance from given S is lower " << endl <<"than minimum possible n*Variance from given X values " << endl;
        //            return -1; // NO FEASIBLE SOLUTION
        }

        // Solve for q to get nVarUbCalc for each loop
        double qtemp1 = 0;
        double qtemp2 = 0;
        double qMeanUb = SpreadConsU.getDomVarUb(0) * SpreadConsX.getDomN();
        double qMeanLb = SpreadConsU.getDomVarLb(0) * SpreadConsX.getDomN();
        // TODO:
        // NOTE: I am not sure if the m = 0 case causes any problems here.
        // Cause if m = 0, the graph can no longer be proven to be convex as the 2nd derivative
        // will be < 0 and is concave instead.

        double a = 0;
        double b = 0; // note: This is b' which is b/2
        double c = 0;
        double inRoot = 0;
        for (int i = 0; i < SpreadConsX.numI; i++)
        {
            // Note: Always write 1.0 instead of 1 when working with doubles
            // if not it will divide as integers and result in an integer answer
            a =((1.0/SpreadConsX.I[i]->m) - (1.0/SpreadConsX.getDomN()));
            b = ((-1.0) * (SpreadConsX.I[i]->ES / SpreadConsX.I[i]->m));
            c = (((SpreadConsX.I[i]->ES * SpreadConsX.I[i]->ES)/((double)SpreadConsX.I[i]->m))+ SpreadConsX.I[i]->C - nVarUbCalc);
            inRoot = b*b - a*c;
            qtemp1 =  ((((-1.0) * b) + sqrt(inRoot))/a);
            qtemp2 = ((((-1.0) * b) - sqrt(inRoot))/a);
// cout << "a: "<< a << " b: " << b << " c: " << c << " inRoot: " << inRoot <<endl;
// cout << "qtemp: " << qtemp1 << "  " << qtemp2 << endl;

            // if the calculated qtemp1 is within the bounds of this interval
            if (qtemp1 >= SpreadConsX.I[i]->Vlb && qtemp1 <= SpreadConsX.I[i]->Vub )
            {
                if (ceil(qtemp1/SpreadConsX.getDomN()) < SpreadConsU.getDomVarUb(0))
                {
                    qMeanUb = qtemp1; // update qMeanUb
                    SpreadConsU.setDomVarUb(0, ceil(qtemp1/SpreadConsX.getDomN()));
                    changed = true;
                    cout << "Success: Pruned U's upper bound to " <<ceil(qtemp1/SpreadConsX.getDomN()) << endl;
                }
            }
            if (qtemp2 >= SpreadConsX.I[i]->Vlb && qtemp2 <= SpreadConsX.I[i]->Vub )
            {
                if (floor(qtemp2/SpreadConsX.getDomN()) > SpreadConsU.getDomVarLb(0))
                {
                    qMeanLb = qtemp2; // update qMeanLb
                    SpreadConsU.setDomVarLb(0, floor(qtemp2/SpreadConsX.getDomN()));
                    changed = true;
                    cout << "Success: Pruned U's lower bound to " <<floor(qtemp2/SpreadConsX.getDomN()) << endl;
                }
            }
        }
        //------------------------------------------------------------

        // Note: This helps in dealing with weird arithmetic error and calculations
        // However, it may result in a less tightly pruned algorithm.
        // It shouldn't matter much in terms of pruning though, cause this case is rare and is often
        // close to ideality.
        if (SpreadConsS.getDomVarUb(0) == 0)
        {
            break; // get out from propagation algorithm for X  if upper bound is 0
        }

        // Prune the X Bounds using S (U is used to determine qSelected)
        // Note: U helps by selecting qSelected
        double qSelected  = 0;
        int indexQSelected = 0;
        int originalIndexQSelected = 0; // to restore original indexQSelected at end of iteration
        int xIndex = 0;
        double DMax = 0;
        int originalLowerBound = 0;
        int maxShift = 0;
        double v = 0;
        double DMin = 0;
        int originalUpperBound = 0;

        // Need these 2 to handle deadling with changing intervals
        // and qSelected
        int XLB[SpreadConsX.getDomN()]; // get all the updated lower bounds
        int XUB[SpreadConsX.getDomN()]; // get all the updated upper bounds

        // For Given Mean Value, when lower bound is equal to upper bound
        if (SpreadConsU.getDomVarLb(0) == SpreadConsU.getDomVarUb(0))
 //if ( 2 == 1) // for debugging else case
//if (1 == 1)
        {
            // note: It is called qSelected here NOT q0, as q0 refers to something else to propagate U below
            qSelected = (double) SpreadConsU.getDomVarLb(0) * SpreadConsX.getDomN();
            //   qSelected = 10; // For Debugging
            cout << "Pruning X using GIVEN U with qSelected: " << qSelected << endl;
            // Determine which interval qSelected is in
            indexQSelected = 0;
            for (int i = 0; i < SpreadConsX.numI; i++)
            {
                if ((qSelected >= SpreadConsX.I[i]->Vlb) && (qSelected <= SpreadConsX.I[i]->Vub))
                {
                    indexQSelected = i;
                    break;
                }
            }
            originalIndexQSelected = indexQSelected;

            xIndex = 0;
            DMax = 0;
            originalLowerBound = 0;
            maxShift = 0;
            // Prune the X variables on R(I) upper bound
            for (int i = 0; i < SpreadConsX.I[originalIndexQSelected]->r; i++)
            {
                xIndex = SpreadConsX.I[originalIndexQSelected]->R[i];
                // Save the originalLowerBound of Current X as it will be changed
                originalLowerBound = SpreadConsX.getDomVarLb(xIndex);
                DMax =  SpreadConsX.findDMax(qSelected,originalIndexQSelected, xIndex, SpreadConsS.getDomVarUb(0));
                cout << "Dmax R(Iq) for X[" << xIndex <<"] is: " << DMax << endl;
                SpreadConsX.setDomVarLb(xIndex, originalLowerBound); // reset original lower bound
                maxShift = ceil(SpreadConsX.getDomVarLb(xIndex) + DMax);
                if (maxShift < SpreadConsX.getDomVarUb(xIndex))
                {
                        // SpreadConsX.setDomVarUb(xIndex, maxShift);
                    XUB[xIndex] = maxShift;
                    changed = true;
                    cout << "Success: Pruned X[" << xIndex << "]'s upper bound to " <<maxShift << endl;
                }
                else
                {
                    XUB[xIndex] = SpreadConsX.getDomVarUb(xIndex);
                }
                // If change upper bound, will not change lower bound,
                XLB[xIndex] = SpreadConsX.getDomVarLb(xIndex);
                // Recreate the intervals to original
                SpreadConsX.createInterval();
            }
            v = 0;
            // Prune the X variables on M(I) upper bound
            // Note: If you did not have orignalIndexQSelected and just indexQSelected,
            // the bottom code will change the value of indexQSelected

            for (int i = 0; i < SpreadConsX.I[originalIndexQSelected]->m ; i++)
            {
                xIndex = SpreadConsX.I[originalIndexQSelected]->M[i];
                cout<< " m is " << xIndex << endl;
                // Save the originalLowerBound of Current X as it will change
                originalLowerBound = SpreadConsX.getDomVarLb(xIndex);
                // Since you can shift all the min to at least up to V, shift it
                // First , get the indexQSelected without shifting
                for (int j = 0; j < SpreadConsX.numI; j++)
                {
                    if ((qSelected >= SpreadConsX.I[j]->Vlb) && (qSelected <= SpreadConsX.I[j]->Vub))
                    {
                        indexQSelected = j; // get new indexQSelected
                        break;
                    }
                }
                // Calculate for v
                v = (qSelected - SpreadConsX.I[indexQSelected]->ES)/((double) SpreadConsX.I[indexQSelected]->m);
                DMax = v - SpreadConsX.getDomVarLb(xIndex); // can shift by at least v - current lower bound
                // Shift current xIndex to v
                SpreadConsX.setDomVarLb(xIndex, floor(v));
                SpreadConsX.createInterval(); // create the new intervals
                // Find the new intervals for indexQSelected
                for (int j = 0; j < SpreadConsX.numI; j++)
                {
                    if ((qSelected >= SpreadConsX.I[j]->Vlb) && (qSelected <= SpreadConsX.I[j]->Vub))
                    {
                        indexQSelected = j;
                        break;
                    }
                }
                DMax +=  SpreadConsX.findDMax(qSelected,indexQSelected, xIndex, SpreadConsS.getDomVarUb(0));
                cout << "Dmax M(Iq) for X[" << xIndex <<"] is: " << DMax << endl;
                SpreadConsX.setDomVarLb(xIndex, originalLowerBound); // reset original lower bound
                maxShift = ceil(SpreadConsX.getDomVarLb(xIndex) + DMax);

                if (maxShift < SpreadConsX.getDomVarUb(xIndex))
                {
                    //SpreadConsX.setDomVarUb(xIndex, maxShift);
                    XUB[xIndex] = maxShift;
                    changed = true;
                    cout << "Success: Pruned X[" << xIndex << "]'s upper bound to " <<maxShift << endl;
                }
                else
                {
                    XUB[xIndex] = SpreadConsX.getDomVarUb(xIndex);
                }
                // If X is in M(Iq), will need to prune both ways
                // Reset the  intervals to original
                SpreadConsX.createInterval();
            }

            DMin = 0;
            originalUpperBound = 0;
            // Prune the X variables on L(I) lower bound
            for (int i = 0; i < SpreadConsX.I[originalIndexQSelected]->l; i++)
            {
                xIndex = SpreadConsX.I[originalIndexQSelected]->L[i];
                originalUpperBound = SpreadConsX.getDomVarUb(xIndex);
                DMin =  SpreadConsX.findDMin(qSelected,originalIndexQSelected, xIndex, SpreadConsS.getDomVarUb(0));
                cout << "Dmin L(Iq) for X[" << xIndex <<"] is: " << DMax << endl;
                SpreadConsX.setDomVarUb(xIndex, originalUpperBound); // reset original lower bound
                maxShift = floor(SpreadConsX.getDomVarUb(xIndex) - DMin);
                if (maxShift > SpreadConsX.getDomVarLb(xIndex))
                {
                    XLB[xIndex] = maxShift;
                    //SpreadConsX.setDomVarLb(xIndex, maxShift);
                    changed = true;
                    cout << "Success: Pruned X[" << xIndex << "]'s lower bound to " <<maxShift << endl;
                }
                else
                {
                    XLB[xIndex] =  SpreadConsX.getDomVarLb(xIndex);
                }
                XUB[xIndex] = SpreadConsX.getDomVarUb(xIndex);
                // Reset to original intervals
                SpreadConsX.createInterval();
            }// end of L algorithm

            v = 0;
            DMin = 0;
            originalUpperBound = 0;
            // Prune the X variables on M(I) lower bound
            for (int i = 0; i < SpreadConsX.I[originalIndexQSelected]->m ; i++)
            {
                xIndex = SpreadConsX.I[originalIndexQSelected]->M[i];
                cout<< " m is " << xIndex << endl;
                // Save the originalLowerBound of Current X as it will change
                originalUpperBound = SpreadConsX.getDomVarUb(xIndex);
                // Since you can shift all the min to at least up to V, shift it
                // First , get the indexQSelected without shifting
                // Calculate for v
                v = (qSelected - SpreadConsX.I[originalIndexQSelected]->ES)/((double) SpreadConsX.I[originalIndexQSelected]->m);
                DMin = SpreadConsX.getDomVarUb(xIndex) - v; // can shift by at least current upper bound - v to v
                // Shift current xIndex to v
                SpreadConsX.setDomVarUb(xIndex, ceil(v));
                SpreadConsX.createInterval(); // create the new intervals
                // Find the new intervals for indexQSelected
                for (int j = 0; j < SpreadConsX.numI; j++)
                {
                    if ((qSelected >= SpreadConsX.I[j]->Vlb) && (qSelected <= SpreadConsX.I[j]->Vub))
                    {
                        indexQSelected = j;
                        break;
                    }
                }
                DMin +=  SpreadConsX.findDMin(qSelected,indexQSelected, xIndex, SpreadConsS.getDomVarUb(0));
                cout << "Dmin M(Iq) for X[" << xIndex <<"] is: " << DMin << endl;
                SpreadConsX.setDomVarUb(xIndex, originalUpperBound); // reset original upper bound
                maxShift = floor(SpreadConsX.getDomVarUb(xIndex) - DMin);

                if (maxShift  > SpreadConsX.getDomVarLb(xIndex))
                {
                    //SpreadConsX.setDomVarUb(xIndex, maxShift);
                    XLB[xIndex] = maxShift;
                    changed = true;
                    cout << "Success: Pruned X[" << xIndex << "]'s lower bound to " <<maxShift << endl;
                }
                else
                {
                    XLB[xIndex] = SpreadConsX.getDomVarLb(xIndex);
                }
                // If X is in M(Iq), will need to prune both ways
                // Reset the  intervals to original
                SpreadConsX.createInterval();
            }

            // Note: Important TODO:
            // Doing this below that you only update values after checking all X
            // is incorrect as you may over-prune but it makes the code a lot faster than
            // simply re-calculating everthing as soon as you find one single X that gets to be pruned.

            // Update all bound values
            for(i = 0; i < SpreadConsX.getDomN(); i++)
            {
                SpreadConsX.setDomVarLb(i, XLB[i]);
                SpreadConsX.setDomVarUb(i, XUB[i]);
            }
            // Recreate the new intervals
            SpreadConsX.createInterval();
            // Re-print latest bounds
            cout<<"Latest Bounds:" << endl;
            SpreadConsX.printDomVar();
            SpreadConsU.printDomVar();
            SpreadConsS.printDomVar();
        } // End of given mean algorithm

        // For general mean values, when lower bound not equal to upper bound of mean
        // Note: The paper is wrong as the DMax(q) is shown to be convex and not derivable from Matlab Code
        // I have confirmed this error with the author himself via email.
        else
        {
            cout << "Pruning X using GENERAL U" << endl;
            qSelected  = SpreadConsU.getDomVarLb(0); // Initialize to the lower bound4
            qSelected = 0; // to debug below
            indexQSelected = 0;
            originalIndexQSelected = 0; // to restore original indexQSelected at end of iteration
            xIndex = 0;
            DMax = 0;
            originalLowerBound = 0;
            maxShift = 0;
            v = 0;
            DMin = 0;
            originalUpperBound = 0;
            double XMin = 0; // current XMin used to calculate DMax for L(Iq), R(Iq)
            double XMax = 0; // current XMax used to calculate DMin for L(Iq), R(Iq)
            double QMax = 0; // Qmax from current U's upper bound
            double QMin = 0; // Qmin from current U's lower bound
            double stdDevUb = SpreadConsS.getDomVarUb(0);
            double optMax = 1.0*stdDevUb*stdDevUb*SpreadConsX.getDomN();
            QMax = SpreadConsU.getDomVarUb(0) * SpreadConsX.getDomN();
            QMin = SpreadConsU.getDomVarLb(0) * SpreadConsX.getDomN();
            stdDevUb = SpreadConsS.getDomVarUb(0);
            double dMaxTemp = 0; // to get value of dMax from DMAX() function
            double dMaxFinal = 0; // initialize to minimum possible value
            double dMinTemp = 0;  // to get value of dMin from DMIN() function
            double dMinFinal = 0; // initialize to minimum possible value
            double qArgumentLB = 0;
            double qArgumentUB = 0;
            bool pruneUseR = false;
            bool pruneUseL = false;
            bool pruneUseM = false;
            bool stopPruning = false;
            // Loop through each variable and calculate qMax to be used fpr R, L, M or none
            for(int s = 0; s < SpreadConsX.getDomN(); s++)
            {
                // Note: Each X can either belong to (R or M) & (L or M) but not both
                // Stop pruning is used to not bother looking for X in M if already prune for R for DMax case
                pruneUseR = false;
                pruneUseM = false;
                stopPruning = false;
                pruneUseL = false;

                XMin = SpreadConsX.getDomVarLb(s);
                XMax = SpreadConsX.getDomVarUb(s);
                xIndex = s;

                //--------------------------------------------------------------------------------------
                // Prune Upper Bound of X using R(Iq) and M(Iq)
                //--------------------------------------------------------------------------------------
                // Loop through every interval and keep track of which interval has the largest dMaxTemp
                for(int t = 0; t < SpreadConsX.numI; t++)
                {
                    qArgumentLB = SpreadConsX.I[t]->Vlb;
                    qArgumentUB = SpreadConsX.I[t]->Vub;
                    // Case 1: Refer to left end point (lower bound)
                    // Note: If you sketch the graph, you can have 6 different scenarios,
                    // each scenario below covers 2 of the 6 by these inequality tests
                    // Scenario 1:  (qArgumentLB >= QMin) && (qArgumentLB <= QMax)
                                    // set qArgumentLB to qArgumentLB
                    // Scenario 2: ( qArgumentLB<= QMin) && (qArgumentUB >= QMin)
                                    // set qArgumentLB to QMin
                    // Scenario 3:  (qArgumentUB < QMin) || (qArgumentLB > QMax)
                                    // continue; cause nothing can be done (including for Case 2)
                    // Scenario 1
                    if((qArgumentLB >= QMin) && (qArgumentLB <= QMax))
                    {
                        // Do nothing cause already assigned correctly
                    }
                    // Scenario 2
                    else if((qArgumentLB<= QMin) && (qArgumentUB >= QMin))
                    {
                        qArgumentLB = QMin;
                    }
                    // Scenario 3
                    else if ((qArgumentUB < QMin) || (qArgumentLB > QMax))
                    {
                        continue; // go to next iteration of loop
                    }
                    else
                    {
                        cout << " SHOULDN'T BE IN HERE" << endl;
                        continue;
                    }
                    dMaxTemp = SpreadConsX.DMAX(qArgumentLB,XMin,optMax,0);
                    cout << "dMaxTemp: " << dMaxTemp << " dMaxFinal: " << dMaxFinal << endl;
                    if (dMaxTemp > dMaxFinal)
                    {
                        dMaxFinal = dMaxTemp;
                        qSelected = qArgumentLB;
                    }
                    // Case 2: Refer to right end point (upper bound)
                    // Scenario 1:  (qArgumentUB >= QMin) && (qArgumentUB <= QMax)
                                    // set qArgumentUB to qArgumentUB
                    // Scenario 2: (qArgumentUB >= QMax) && (qArgumentLB <= QMax)
                                    // set qArgumentLB to QMin
                    // Scenario 3:  (qArgumentUB < QMin) || (qArgumentLB > QMax)

                    // Scenario 1
                    if ((qArgumentUB >= QMin) && (qArgumentUB <= QMax))
                    {
                        // do nothing, qArgumentUB is already the right value
                    }
                    // Scenario 2
                    else if((qArgumentUB >= QMax) && (qArgumentLB <= QMax))
                    {
                        qArgumentUB = QMax;
                    }
                    // Scenario 3
                    else if ((qArgumentUB < QMin) || (qArgumentLB > QMax))
                    {
                        cout << " Should never be in here since above Case 1 should have continue to next iteration" << endl;
                        continue;
                    }
                    else
                    {
                        cout << " Definitely should never be in here" << endl;
                        continue;
                    }
                    dMaxTemp = SpreadConsX.DMAX(qArgumentUB,XMin,optMax,1);
                    if (dMaxTemp > dMaxFinal)
                    {
                        dMaxFinal = dMaxTemp;
                        qSelected = qArgumentUB;
                    }
                } // End of current interval
                // Here, have finish assigning dMaxFinal and qSelected for R(Iq) and M(Iq)
                // To prune upper bound of X
                // Need calculate dMax from qSelected to prune R(Iq) and M(Iq)
                // First, get QIndex
                cout << "Using DMAX to prune upper bound of current X" << endl;
                cout << "xIndex is " << xIndex << endl;
                cout << "qSelected is " << qSelected << endl;
                // Since you can shift all the min to at least up to V, shift it
                // First , get the indexQSelected without shifting
                for (int j = 0; j < SpreadConsX.numI; j++)
                {
                    if ((qSelected >= SpreadConsX.I[j]->Vlb) && (qSelected <= SpreadConsX.I[j]->Vub))
                    {
                        indexQSelected = j; // get new indexQSelected
                        break;
                    }
                }
                cout << "indexQSelected is " << indexQSelected << endl;
                originalIndexQSelected = indexQSelected;
                for (int i = 0; i < SpreadConsX.I[originalIndexQSelected]->r; i++)
                {
                    // If this X is part of R
                    if (xIndex == SpreadConsX.I[originalIndexQSelected]->R[i])
                    {
                        pruneUseR = true;
                        stopPruning = true; // don't bother pruning M for this X
                    }
                }
                // Step 1: Prune this X using R(Iq) if X is in R(Iq)
                if (pruneUseR)
                {
                    cout << "Pruning upperbound using R(Iq)" << endl;
                    // Save the originalLowerBound of Current X as it will change
                    originalIndexQSelected = indexQSelected;
                    DMax = 0;
                    originalLowerBound = 0;
                    maxShift = 0;
                    // Prune current X on R(I) upper bound
                    // Save the originalLowerBound of Current X as it will be changed
                    originalLowerBound = SpreadConsX.getDomVarLb(xIndex);
                    DMax =  SpreadConsX.findDMax(qSelected,originalIndexQSelected, xIndex, SpreadConsS.getDomVarUb(0));
                    cout << "Dmax R(Iq) for X[" << xIndex <<"] is: " << DMax << endl;
                    SpreadConsX.setDomVarLb(xIndex, originalLowerBound); // reset original lower bound
                    maxShift = ceil(SpreadConsX.getDomVarLb(xIndex) + DMax);
                    if (maxShift < SpreadConsX.getDomVarUb(xIndex))
                    {
                        // SpreadConsX.setDomVarUb(xIndex, maxShift);
                        XUB[xIndex] = maxShift;
                        changed = true;
                        cout << "Success: Pruned X[" << xIndex << "]'s upper bound to " <<maxShift << endl;
                    }
                    else
                    {
                        XUB[xIndex] = SpreadConsX.getDomVarUb(xIndex);
                    }
                    // Recreate the intervals to original
                    SpreadConsX.createInterval();
                } // end of if(pruneUseR)

                // If X wasn't in R
                if(stopPruning == false)
                {
                    for (int i = 0; i < SpreadConsX.I[originalIndexQSelected]->m; i++)
                    {
                        // If this X is part of M
                        if (xIndex == SpreadConsX.I[originalIndexQSelected]->M[i])
                        {
                            pruneUseM = true;
                            stopPruning = true; // don't setting default values for XUB for this X
                        }
                    }
                }
                // Step 2: Prune this X using M(Iq)
                if (pruneUseM)
                {
                    cout << "Pruning upperbound using M(Iq)" << endl;

                    v = 0;
                    // Save the originalLowerBound of Current X as it will change
                    originalLowerBound = SpreadConsX.getDomVarLb(xIndex);
                    // Since you can shift all the min to at least up to V, shift it
                    // First , get the indexQSelected without shifting
                    // If in here, m won't be 0, so can definitely calculate for v
                    // Calculate for v
                    v = (qSelected - SpreadConsX.I[originalIndexQSelected]->ES)/((double) SpreadConsX.I[originalIndexQSelected]->m);
                    DMax = v - SpreadConsX.getDomVarLb(xIndex); // can shift by at least v - current lower bound
                    // Shift current xIndex to v
                    SpreadConsX.setDomVarLb(xIndex, floor(v));
                    SpreadConsX.createInterval(); // create the new intervals
                    // Find the new intervals for indexQSelected
                    for (int i = 0; i< SpreadConsX.numI; i++)
                    {
                        if ((qSelected >= SpreadConsX.I[i]->Vlb) && (qSelected <= SpreadConsX.I[i]->Vub))
                        {
                            indexQSelected = i;
                            break;
                        }
                    }
                    DMax +=  SpreadConsX.findDMax(qSelected,indexQSelected, xIndex, SpreadConsS.getDomVarUb(0));
                    cout << "Dmax M(Iq) for X[" << xIndex <<"] is: " << DMax << endl;
                    SpreadConsX.setDomVarLb(xIndex, originalLowerBound); // reset original lower bound
                    maxShift = ceil(SpreadConsX.getDomVarLb(xIndex) + DMax);

                    if (maxShift < SpreadConsX.getDomVarUb(xIndex))
                    {
                        //SpreadConsX.setDomVarUb(xIndex, maxShift);
                        XUB[xIndex] = maxShift;
                        changed = true;
                        cout << "Success: Pruned X[" << xIndex << "]'s upper bound to " <<maxShift << endl;
                    }
                    else
                    {
                        XUB[xIndex] = SpreadConsX.getDomVarUb(xIndex);
                    }
                    // Reset the  intervals to original
                    SpreadConsX.createInterval();
                }// end of if(pruneUseM)

                // If did not prune using R or M for this X, set to original value
                if(stopPruning == false)
                {
                    XUB[xIndex] = SpreadConsX.getDomVarUb(xIndex);
                }

                //--------------------------------------------------------------------------------------
                // Prune Lower Bound of X using L(Iq) and M(Iq)
                //--------------------------------------------------------------------------------------
                // Find qSelect and indexQSelect using DMIN
                // Loop through every interval and keep track of which interval has the largest dMaxTemp

                // Note: Each X can either belong to (R or M) & (L or M) but not both
                // Stop pruning is used to not bother looking for X in M if already prune for R for DMax case
                pruneUseL = false;
                pruneUseM = false;
                stopPruning = false;

                for(int t = 0; t < SpreadConsX.numI; t++)
                {
                    qArgumentLB = SpreadConsX.I[t]->Vlb;
                    qArgumentUB = SpreadConsX.I[t]->Vub;
                    // Case 1: Refer to left end point (lower bound)
                    // Note: If you sketch the graph, you can have 6 different scenarios,
                    // each scenario below covers 2 of the 6 by these inequality tests
                    // Scenario 1:  (qArgumentLB >= QMin) && (qArgumentLB <= QMax)
                                    // set qArgumentLB to qArgumentLB
                    // Scenario 2: ( qArgumentLB<= QMin) && (qArgumentUB >= QMin)
                                    // set qArgumentLB to QMin
                    // Scenario 3:  (qArgumentUB < QMin) || (qArgumentLB > QMax)
                                    // continue; cause nothing can be done (including for Case 2)
                    // Scenario 1
                    if((qArgumentLB >= QMin) && (qArgumentLB <= QMax))
                    {
                        // Do nothing cause already assigned correctly
                    }
                    // Scenario 2
                    else if((qArgumentLB<= QMin) && (qArgumentUB >= QMin))
                    {
                        qArgumentLB = QMin;
                    }
                    // Scenario 3
                    else if ((qArgumentUB < QMin) || (qArgumentLB > QMax))
                    {
                        continue; // go to next iteration of lSpreadConsXoop
                    }
                    else
                    {
                        cout << " SHOULDN'T BE IN HERE" << endl;
                        continue;
                    }
                    dMinTemp = SpreadConsX.DMIN(qArgumentLB,XMax,optMax,0);
                    if (dMinTemp > dMinFinal)
                    {
                        dMinFinal = dMinTemp;
                        qSelected = qArgumentLB;
                    }
                    // Case 2: Refer to right end point (upper bound)
                    // Scenario 1:  (qArgumentUB >= QMin) && (qArgumentUB <= QMax)
                                    // set qArgumentUB to qArgumentUB
                    // Scenario 2: (qArgumentUB >= QMax) && (qArgumentLB <= QMax)
                                    // set qArgumentLB to QMin
                    // Scenario 3:  (qArgumentUB < QMin) || (qArgumentLB > QMax)

                    // Scenario 1
                    if ((qArgumentUB >= QMin) && (qArgumentUB <= QMax))
                    {
                        // do nothing, qArgumentUB is already the right value
                    }
                    // Scenario 2
                    else if((qArgumentUB >= QMax) && (qArgumentLB <= QMax))
                    {
                        qArgumentUB = QMax;
                    }
                    // Scenario 3
                    else if ((qArgumentUB < QMin) || (qArgumentLB > QMax))
                    {
                        cout << " Should never be in here since above Case 1 should have continue to next iteration" << endl;
                        continue;
                    }
                    else
                    {
                        cout << " Definitely should never be in here" << endl;
                        continue;
                    }
                    dMinTemp = SpreadConsX.DMIN(qArgumentUB,XMax,optMax,1);
                    if (dMinTemp > dMinFinal)
                    {
                        dMinFinal = dMinTemp;
                        qSelected = qArgumentUB;
                    }
                } // End of current interval
                // Here, have finish assigning dMinFinal and qSelected for L(Iq) and M(Iq)
                // To prune lower bound of X
                // Need calculate dMin from qSelected to prune L(Iq) and M(Iq)
                // First, get QIndex
                cout << "Using DMin to prune lower bound of current X" << endl;
                cout << "xIndex is " << xIndex << endl;
                cout << "qSelected is " << qSelected << endl;
                // Since you can shift all the min to at least up to V, shift it
                // First , get the indexQSelected without shifting
                for (int j = 0; j < SpreadConsX.numI; j++)
                {
                    if ((qSelected >= SpreadConsX.I[j]->Vlb) && (qSelected <= SpreadConsX.I[j]->Vub))
                    {
                        indexQSelected = j; // get new indexQSelected
                        break;
                    }
                }
                cout << "indexQSelected is " << indexQSelected << endl;
                originalIndexQSelected = indexQSelected;
                for (int i = 0; i < SpreadConsX.I[originalIndexQSelected]->l; i++)
                {
                    // If this X is part of L
                    if (xIndex == SpreadConsX.I[originalIndexQSelected]->L[i])
                    {
                        pruneUseL = true;
                        stopPruning = true; // don't bother pruning M for this X
                    }
                }
                // Step 3: Prune this X using L(Iq) if X is in L(Iq)
                if (pruneUseL)
                {
                    cout << "Pruning lowerbound using L(Iq)" << endl;
                    // Save the originalLowerBound of Current X as it will change
                    DMin = 0;
                    originalUpperBound = 0;
                    maxShift = 0;
                    // Prune current X on L(I) upper bound
                    // Save the originalUpperBound of Current X as it will be changed
                    originalUpperBound = SpreadConsX.getDomVarUb(xIndex);
                    DMin =  SpreadConsX.findDMin(qSelected,originalIndexQSelected, xIndex, SpreadConsS.getDomVarUb(0));
                    cout << "Dmin L(Iq) for X[" << xIndex <<"] is: " << DMin << endl;
                    SpreadConsX.setDomVarUb(xIndex, originalUpperBound); // reset original upper bound
                    maxShift = floor(SpreadConsX.getDomVarUb(xIndex) - DMin);
                    if (maxShift > SpreadConsX.getDomVarLb(xIndex))
                    {
                        // SpreadConsX.setDomVarLb(xIndex, maxShift);
                        XLB[xIndex] = maxShift;
                        changed = true;
                        cout << "Success: Pruned X[" << xIndex << "]'s lower bound to " << maxShift << endl;
                    }
                    else
                    {
                        XLB[xIndex] = SpreadConsX.getDomVarLb(xIndex);
                    }
                    // Recreate the intervals to original
                    SpreadConsX.createInterval();
                } // end of if(pruneUseR)

                // If X wasn't in R
                if(stopPruning == false)
                {
                    for (int i = 0; i < SpreadConsX.I[originalIndexQSelected]->m; i++)
                    {
                        // If this X is part of M
                        if (xIndex == SpreadConsX.I[originalIndexQSelected]->M[i])
                        {
                            pruneUseM = true;
                            stopPruning = true; // don't setting default values for XUB for this X
                        }
                    }
                }
                // Step 4: Prune this X using M(Iq) if X is in M(Iq)
                if (pruneUseM)
                {
                    cout << "Pruning lowerbound using M(Iq)" << endl;
                    v = 0;
                    // Save the originalUpperBound of Current X as it will change
                    originalUpperBound = SpreadConsX.getDomVarUb(xIndex);
                    // Since you can shift all the max to at least down to V, shift it
                    // First , get the indexQSelected without shifting
                    // Note: If in here, definitely m isn't 0, so can calculate for v normally
                    // Calculate for v
                    v = (qSelected - SpreadConsX.I[originalIndexQSelected]->ES)/((double) SpreadConsX.I[originalIndexQSelected]->m);
                    DMin = SpreadConsX.getDomVarUb(xIndex) - v; // can shift by at least current upper bound - v to v
                    // Shift current xIndex to v
                    SpreadConsX.setDomVarUb(xIndex, ceil(v));
                    SpreadConsX.createInterval(); // create the new intervals
                    // Find the new intervals for indexQSelected
                    for (int i = 0; i < SpreadConsX.numI; i++)
                    {
                        if ((qSelected >= SpreadConsX.I[i]->Vlb) && (qSelected <= SpreadConsX.I[i]->Vub))
                        {
                            indexQSelected = i;
                            break;
                        }
                    }
                    DMin +=  SpreadConsX.findDMin(qSelected,indexQSelected, xIndex, SpreadConsS.getDomVarUb(0));
                    cout << "Dmin M(Iq) for X[" << xIndex <<"] is: " << DMin << endl;
                    SpreadConsX.setDomVarUb(xIndex, originalUpperBound); // reset original upper bound
                    maxShift = floor(SpreadConsX.getDomVarUb(xIndex) - DMin);
                    if (maxShift  > SpreadConsX.getDomVarLb(xIndex))
                    {
                        //SpreadConsX.setDomVarUb(xIndex, maxShift);
                        XLB[xIndex] = maxShift;
                        changed = true;
                        cout << "Success: Pruned X[" << xIndex << "]'s lower bound to " <<maxShift << endl;
                    }
                    else
                    {
                        XLB[xIndex] = SpreadConsX.getDomVarLb(xIndex);
                    }
                    // If X is in M(Iq), will need to prune both ways
                    // Reset the  intervals to original
                    SpreadConsX.createInterval();
                }// end of if(pruneUseM)

                // If did not prune using L or M for this X, set to original value
                if(stopPruning == false)
                {
                    XLB[xIndex] = SpreadConsX.getDomVarLb(xIndex);
                }
            } // End of current X

            // Note: Important : TODO:
            // Doing this below that you only update values after checking all X
            // is incorrect as you may over-prune but it makes the code a lot faster than
            // simply re-calculating everthing as soon as you find one single X that gets to be pruned.

            // After done checking all X
            // Update all bound values
            for(i = 0; i < SpreadConsX.getDomN(); i++)
            {
                SpreadConsX.setDomVarLb(i, XLB[i]);
                SpreadConsX.setDomVarUb(i, XUB[i]);
            }
            // Recreate the new intervals
            SpreadConsX.createInterval();
            // Re-print latest bounds
            cout<<"Latest Bounds:" << endl;
            SpreadConsX.printDomVar();
            SpreadConsU.printDomVar();
            SpreadConsS.printDomVar();
        } // End of Pruning X by General Case of Mean
        //------------------------------------------------------------
    }
    cout << " SOON SOON SOON SOON" << endl;
    // End of while loop, here, means can't propagate anymore cause nothing is changing
    SpreadConsX.printDomVar();
    SpreadConsS.printDomVar();
    SpreadConsU.printDomVar();
    // End of Propagation algorithm
    int N2 = (consdata->nvars-2); // number of Variables for X 

    // In SCIP , just copy paste everything, and change the main() function to the CONSPROP() propagation function
	for(int i = 0; i < N2; i++)
	{
		SCIP_CALL ( SCIPtightenVarLb(scip, consdata->vars[i], SpreadConsX.getDomVarLb(i), TRUE, &infeasible, NULL ) );
		SCIP_CALL ( SCIPtightenVarUb(scip, consdata->vars[i], SpreadConsX.getDomVarUb(i),, TRUE, &infeasible, NULL ) );  
	}
		SCIP_CALL ( SCIPtightenVarLb(scip, consdata->vars[consdata->nvars-2], SpreadConsU.getDomVarLb(0), TRUE, &infeasible, NULL ) );
		SCIP_CALL ( SCIPtightenVarUb(scip, consdata->vars[consdata->nvars-2], SpreadConsU.getDomVarUb(0), TRUE, &infeasible, NULL ) );
		SCIP_CALL ( SCIPtightenVarLb(scip, consdata->vars[consdata->nvars-1], SpreadConsS.getDomVarLb(0), TRUE, &infeasible, NULL ) );
		SCIP_CALL ( SCIPtightenVarUb(scip, consdata->vars[consdata->nvars-1], SpreadConsS.getDomVarUb(0), TRUE, &infeasible, NULL ) );

		// DONE COPYING PASTING ANSWER BACK INTO SCIP 
   
   /* THIAGO's Implementation 
     	SCIP_CONSDATA* consdata;
	SCIP_CONS* cons;
     	cons = conss[c];
     	assert( cons != NULL );	
     	consdata = SCIPconsGetData(cons);

	// Get number of variables 
        int nvars = (consdata->nvars-2); // number of Variables for X 

//TODO: means must have bounds and q would not just be given a value from u 
	// get q = number of variables * mean 
// TODO: Now is just always using upper bound, change to getting optimal value here
        double q = SCIPvarGetUbLocal(consdata->vars[consdata->nvars-1])   *nvars;
// TODO: IF statement to check if variable mean or any mean given 
	// Get Sub and Slb  
        double sum_min = 0;
        double sum_max = 0;
 	double vars[nvars][2]; // X variables 
        for (int k = 0; k < nvars; k++)
        {
 	     assert(consdata->vars[k] != NULL);
             vars[k][0] = SCIPvarGetLbLocal(consdata->vars[k]);
             vars[k][1] = SCIPvarGetUbLocal(consdata->vars[k]);
			// printf("\n (OLD) Name = %s, Lower bound = %f, Upper bound = %f.", SCIPvarGetName(consdata->vars[k]), vars[k][0], vars[k][1]);
             sum_min += vars[k][0];
             sum_max += vars[k][1];
        }

	// Q must be an element of [Slb, Sub], if it is not, return SCIP_CUTOFF 
        //checks if q belongs to set of feasible solutions
        if ((sum_min > q) || (sum_max < q))
        {
	   *result = SCIP_CUTOFF;
           return SCIP_OKAY;
        }

    	//calculates an upper bound on the std dev based on domains of the variable
        double ub_stddev = max_stddev(nvars, vars, SCIPvarGetUbLocal(consdata->vars[consdata->nvars-2])); // using helper function defined above (pg 13 of 15) 
	// Note: numerical_gap is a const. defined at the top of this file
	// Question: No idea why need numerical_gap, prolly some numerical methods problem 
	// If the max std. deviation is below the current upper bound, tighten it. 
	if (ub_stddev + numerical_gap < SCIPvarGetUbLocal(consdata->vars[consdata->nvars-1])) // get maxStdDev
	   SCIP_CALL ( SCIPtightenVarUb	(scip, consdata->vars[consdata->nvars-1], ub_stddev, TRUE, &infeasible, NULL ) );

	// Maximum std dev criterion
        double upper_bound = SCIPvarGetUbLocal(consdata->vars[consdata->nvars-1]); // upper bound of std deviation 
        double pi1Max = nvars * pow ( upper_bound, 2.0); //  pi = numVar * variance 


	//set of values B, maximum of (2n - 1) number of intervals (refer to pg 4 of 15) 
        double *bounds = new double[nvars*2+5]; // Question: No idea why he added 5 here 
	int nb;
	create_bounds(nvars, vars, nb, bounds);// Create the contiguous intervals (refer pg 3 of 15) 

	//intervals I
	int nI = nb-1; // number of intervals in I is 1 less than number of bounds in B , e.g. 2,3,4 => 3 bounds, but 2 intervals [2,3], [3,4]
	interval* I = new interval[nI+5]; // Create an array of intervals // Question: No idea why he added 5 here again 
	for (int j=0;j<nI;j++)	  
	{
	   I[j].min = bounds[j];
	   I[j].max = bounds[j+1];	
	   calculate_indices(nvars, vars, I[j]); // Calculate the other properties such as ES, R, L, M, and V (pg 5 of 15) 
	}

// TODO: INSERT UPDATED CODE HERE 


// TODO: NOTE: EVERYTHING BELOW HERE ASSUMES q IS GIVEN, YOU NEED TO CALCULATE A VALUE FOR U BEFORE YOU CAN USE ANYTHING BELOW HERE. 
	//propagator algorithm
	int changes = 0; //0: no changes, 1: changes
        double d;
	double min, max; //store original values of the variables before the shiftings
	for (int j=0;j < nI;j++) 
	{
	   // Loop through all intervals and allow execute if q fits within this interval 
           //if the q fits in the interval V(Iq), we select Iq
	   if ((I[j].V_min <= q) && (q <= I[j].V_max)) 
	   { 
   	      //variables to the Right of the interval j
 	      for (int i=0;i < I[j].r; i++)
              { 
		  //if the lower bound is equal to the upper bound, there's no possible improvement
		  if (vars[I[j].R[i]][0] == vars[I[j].R[i]][1])
			continue; // get out of current Right X variable and move to next Right X variable 
        	  //stores original values of the variables cause FindDMax function changes their bounds 
	          min = vars[I[j].R[i]][0];
        	  max = vars[I[j].R[i]][1];
	          d = FindDMax(nvars, vars, I[j].R[i], j, pi1Max, q, I[j]);
        	  vars[I[j].R[i]][0] = min;
	          vars[I[j].R[i]][1] = max;
		  // Note: If shift is negative, it means problem is infeasible at current node as it has to shift below the minimum bound. 
        	  if (d < 0) 
	          {
		       // Clean up memory handling 
		       delete [] bounds;
                       for (int k=0;k < nI; k++)
 	 	       {
	  		  delete [] I[k].L;
			  delete [] I[k].R;
			  delete [] I[k].M;
	 	       }
		       delete [] I;
                       *result = SCIP_CUTOFF;
                       return SCIP_OKAY;
        	  }
	          //If upper bound is more than maximum shift that can occur to the lower bound upwards, tighten the upper bound downwards
  		  if (vars[I[j].R[i]][1] > vars[I[j].R[i]][0] + d + numerical_gap) 
 	          {
		      vars[I[j].R[i]][1] = vars[I[j].R[i]][0] + d;
		      changes = 1;
		  }
              } 
	      //variables to the Left of the interval j       
	      for (int i=0;i < I[j].l; i++)
	      { 
		  //if the lower bound is equal to the upper bound, there's no possible improvement
		  if (vars[I[j].L[i]][0] == vars[I[j].L[i]][1])
			continue;// get out of current Left X variable and move to next Left X variable 
	          //stores original values of the variables
	          min = vars[I[j].L[i]][0];
	          max = vars[I[j].L[i]][1];
	          d = FindDMin(nvars, vars, I[j].L[i], j, pi1Max, q, I[j], nI-1);
	          vars[I[j].L[i]][0] = min;
	          vars[I[j].L[i]][1] = max;
	          if (d < 0) //problem is infeasible
	          {
   	  	       delete [] bounds;
                       for (int k=0;k < nI; k++)
 	 	       {
	  		  delete [] I[k].L;
			  delete [] I[k].R;
			  delete [] I[k].M;
	 	       }
 	               delete [] I;
                       *result = SCIP_CUTOFF;
                       return SCIP_OKAY;
	          }
	          //If lower bound is lower than maximum shift that can occur to the upper bound downwards, tighten the lower bound upwards
	          if (vars[I[j].L[i]][0] < vars[I[j].L[i]][1] - d - numerical_gap) 
		  {
    	             vars[I[j].L[i]][0] = vars[I[j].L[i]][1] - d;
		     changes = 1;
		  }
              } 
              //variables that contain the interval j (M(I))
	      for (int i=0;i < I[j].m; i++)
	      { 		 
		  //if the lower bound is equal to the upper bound, there's no possible improvement
		  if (vars[I[j].M[i]][0] == vars[I[j].M[i]][1])
			continue; // Skip to next M[I] for I 
	          double v;
		  // refer to (pg 10 of 15) 
	          Calculate_optPI(I[j], q, vars, nvars, v); // store result into v 

	          min = vars[I[j].M[i]][0];
		  max = vars[I[j].M[i]][1];

// QUESTION: DONT UNDERSTAND FROM HERE, W.Y don't understand too, understand yourself ==" 
             //enforce bounds consistency on upper bound
	          vars[I[j].M[i]][1] = vars[I[j].M[i]][1] + v - vars[I[j].M[i]][0];
	          vars[I[j].M[i]][0] = v;
		 //creates new problem after adding new variable x'
	         int nb_;
	         double *bounds_ = new double[2*nvars+5];
	         create_bounds(nvars, vars, nb_, bounds_);
	         int nI_ = nb_-1;
	         interval* I_ = new interval[nI_ + 1+5]; // add 1 extra interval (refer to pg 7 of 15) 
	         int count = 0;
	         int j_;
	         for (int k=count;k < nI_; k++)
	         {
	             I_[k].min = bounds_[k];
	             I_[k].max = bounds_[k+1];	
	             if (bounds_[k+1] == v)
	                j_ = k;
        	 }
	         if (bounds_[0] == v) //if there is no interval completely to the left of the lowest value, we create an interval [v v]
	              j_ = 0;
    	         calculate_indices(nvars, vars, I_[j_]);
                 
   	         d = v - min + FindDMax(nvars, vars, I[j].M[i], j_, pi1Max, q, I_[j_]); // pg 7 of 15
 	         vars[I[j].M[i]][0] = min;
                 vars[I[j].M[i]][1] = max;
                 if (d < 0) //problem is infeasible
	         {
      	               delete [] bounds;
                       for (int k=0;k < nI; k++)
 	 	       {
	  		  delete [] I[k].L;
			  delete [] I[k].R;
			  delete [] I[k].M;
	 	       }
 	               delete [] I;
   	  	       delete [] bounds_;
	  	       delete [] I_[j_].L;
		       delete [] I_[j_].R;
 		       delete [] I_[j_].M;
 	               delete [] I_;
                       *result = SCIP_CUTOFF;
                       return SCIP_OKAY;
        	 }
		 // Tighten upper bound if its more than the maximum shift that can occur to lower bound 
        	 if (vars[I[j].M[i]][1] > vars[I[j].M[i]][0] + d + numerical_gap) 
		 {
        	    vars[I[j].M[i]][1] = vars[I[j].M[i]][0] + d;
                    changes = 1;
		 }

       		 //if the lower bound is equal to the upper bound, there's no possible improvement
		 if (vars[I[j].M[i]][0] == vars[I[j].M[i]][1])
		      continue;
	         //enforce bounds consistency on lower bound
	         min = vars[I[j].M[i]][0];
	         max = vars[I[j].M[i]][1];
	         vars[I[j].M[i]][0] = vars[I[j].M[i]][0] + v - vars[I[j].M[i]][1];
	         vars[I[j].M[i]][1] = v;
	         //creates new problem after adding new variable
	         delete [] bounds_;
	  	 delete [] I_[j_].L;
	         delete [] I_[j_].R;
 		 delete [] I_[j_].M;
	         delete [] I_;

	         bounds_ = new double[2*nvars+5];
	         create_bounds(nvars, vars, nb_, bounds_);
	         nI_ = nb_-1;
	         I_ = new interval[nI_ + 1+5];
	         for (int k=0;k < nI_; k++)
	         {
	             I_[k].min = bounds_[k];
	             I_[k].max = bounds_[k+1];	
	             if (bounds_[k] == v)
	                j_ = k;
	         }
	         if (bounds_[nI_] == v) //if there is no interval completely to the right of the highest value, we create an interval [v v]
	            j_ = nI_-1;
                 calculate_indices(nvars, vars, I_[j_]);
   	         d = - v + max + FindDMin(nvars, vars, I[j].M[i], j_, pi1Max, q, I_[j_], nI_-1);
 	         vars[I[j].M[i]][0] = min;
	         vars[I[j].M[i]][1] = max;
	         if (d < 0)  //problem is infeasible
	         {
		       delete [] bounds;
                       for (int k=0;k < nI; k++)
 	 	       {
	  		  delete [] I[k].L;
			  delete [] I[k].R;
			  delete [] I[k].M;
	 	       }
 	               delete [] I;
   	  	       delete [] bounds_;
	  	       delete [] I_[j_].L;
	               delete [] I_[j_].R;
 		       delete [] I_[j_].M;
	               delete [] I_;
                       *result = SCIP_CUTOFF;
                       return SCIP_OKAY;
	         }
		// Tighten lower bound 
        	 if (vars[I[j].M[i]][0] < vars[I[j].M[i]][1] - d - numerical_gap) 
		 {
	            vars[I[j].M[i]][0] = vars[I[j].M[i]][1] - d;
    	            changes = 1;
		 }
	         delete [] bounds_;
	  	 delete [] I_[j_].L;
	         delete [] I_[j_].R;
 		 delete [] I_[j_].M;
                 delete [] I_;
	      } // End of for loop for M(I) 
	      j = nI;
	   } // End of if q is in I 
	} // end of looping through I 

	// Clean up 
        delete [] bounds;
        for (int k=0;k < nI; k++)
        {
           delete [] I[k].L;
           delete [] I[k].R;
	   delete [] I[k].M;
	}
        delete [] I;

	// If changes were made, just loop through all variables and tighten bounds. 
	// Don't really like this implementation cause it has to go through everything although it makes it easier to code. 
	if (changes)
	{
	    printf("Latest bounds are\n");
            for (int k = 0; k < nvars; k++)
            {
//               printf("\n (NEW) Name = %s, Lower bound = %f, Upper bound = %f.", SCIPvarGetName(consdata->vars[k]), vars[k][0], vars[k][1]);
	       SCIP_CALL ( SCIPtightenVarLb	(scip, consdata->vars[k], vars[k][0], TRUE, &infeasible, NULL ) );
  	       SCIP_CALL ( SCIPtightenVarUb	(scip, consdata->vars[k], vars[k][1], TRUE, &infeasible, NULL ) );
		printf("X[ %d ] LB= %f  UB= %f \n", k, SCIPvarGetLbLocal(consdata->vars[k]), SCIPvarGetUbLocal(consdata->vars[k])); 
	    }
	    printf("STDDEV LB= %f  UB= %f \n", SCIPvarGetLbLocal(consdata->vars[nvars]), SCIPvarGetUbLocal(consdata->vars[nvars])); 
	    printf("MEAN LB= %f  UB= %f \n", SCIPvarGetLbLocal(consdata->vars[nvars+1]), SCIPvarGetUbLocal(consdata->vars[nvars+1])); 
	    *result = SCIP_REDUCEDDOM;
     	    return SCIP_OKAY;
        }

		
*/
		
   }

   return SCIP_OKAY;
}



/** presolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRESOL(consPresolSpread)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Spread constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolSpread NULL
#endif


/** propagation conflict resolving method of constraint handler */
//#if 0
static
SCIP_DECL_CONSRESPROP(consRespropSpread)
{  /*lint --e{715}*/
   SCIPdebugMessage("\nConsresprop");
   SCIPerrorMessage("method of Spread constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
//#else
//#define consRespropSpread NULL
//#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockSpread)
{  
    SCIP_CONSDATA* consdata;
    SCIP_VAR** vars;

    SCIPdebugMessage("lock Spread constraint <%s> with nlockspos = %d, nlocksneg = %d\n", SCIPconsGetName(cons), nlockspos, nlocksneg);

    assert(scip != NULL);
    assert(cons != NULL);

    SCIPdebugMessage("\nConsLock");


    consdata = SCIPconsGetData(cons);
    assert(consdata != NULL);

    int i;

    vars = consdata->vars;
    assert(vars != NULL);

    /* In a general way, a variable can make the solution infeasible if rounded both up and down. Here, we tell scip not to round variables in both ways. */
    for( i = 0; i < consdata->nvars; i++) // Note: Include locking mean and std deviation 
    {
        SCIP_CALL( SCIPaddVarLocks(scip, vars[i], nlockspos + nlocksneg, nlockspos + nlocksneg) );
    }
    return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveSpread)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Spread constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}



/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveSpread)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Spread constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** constraint enabling notification method of constraint handler */
static
SCIP_DECL_CONSENABLE(consEnableSpread)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Spread constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** constraint disabling notification method of constraint handler */
static
SCIP_DECL_CONSDISABLE(consDisableSpread)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Spread constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** variable deletion of constraint handler */
static
SCIP_DECL_CONSDELVARS(consDelvarsSpread)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Spread constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/**/
static
SCIP_RETCODE consdataPrint (
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< xor constraint data */
   FILE*                 file               /**< output file (or NULL for standard output) */
   )
{
   assert(consdata != NULL);

/* print coefficients */
   SCIPinfoMessage( scip, file, "Spread({");
   for( int v = 0; v < consdata->nvars; ++v )
   {
      assert(consdata->vars[v] != NULL);
      if( v > 0 )
         SCIPinfoMessage(scip, file, ", ");
      SCIPinfoMessage(scip, file, "%s", SCIPvarGetName(consdata->vars[v]) );
   }
							// mean 								// stdDev 
   SCIPinfoMessage( scip, file, "}, %f, %f])", consdata->vars[consdata->nvars-1], (double)SCIPvarGetUbLocal(consdata->vars[consdata->nvars-2]));
}


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintSpread)
{ 

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   consdataPrint(scip, SCIPconsGetData(cons), file);
   return SCIP_OKAY;

}


/** constraint copying method of constraint handler */
#if 0
static
SCIP_DECL_CONSCOPY(consCopySpread)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Spread constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consCopySpread NULL
#endif


/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseSpread)
{  /*lint --e{715}*/
   SCIPdebugMessage("\nConsParse");
   SCIPerrorMessage("method of Spread constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsSpread)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( varssize < consdata->nvars )
      (*success) = FALSE;
   else
   {
      assert(vars != NULL);
      BMScopyMemoryArray(vars, consdata->vars, consdata->nvars);
      (*success) = TRUE;
   }

   return SCIP_OKAY;}

/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsSpread)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   (*nvars) = consdata->nvars;
   (*success) = TRUE;

   return SCIP_OKAY;

}

/*
 * constraint specific interface methods
 */


/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecSpread)
{  //lint --e{715}
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);

   consdata = (SCIP_CONSDATA*)eventdata;
   assert(consdata != NULL);

   return SCIP_OKAY;
}


/** creates the handler for Spread constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrSpread(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_EVENTHDLR* eventhdlr;

   /* create event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecSpread, NULL) );

   /* create Spread constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata, eventhdlr) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpSpread, consEnfopsSpread, consCheckSpread, consLockSpread,
         conshdlrdata) );


   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopySpread, consCopySpread) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteSpread) );
//#ifdef SCIP_STATISTIC
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreSpread) );
//#endif
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolSpread) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeSpread) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsSpread) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsSpread) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreSpread) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpSpread) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseSpread) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolSpread, CONSHDLR_MAXPREROUNDS, CONSHDLR_DELAYPRESOL) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintSpread) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropSpread, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropSpread) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpSpread, consSepasolSpread, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransSpread) );

   /* add Spread constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/coretimes", "should core-times be propagated (time tabling)?",
         &conshdlrdata->coretimes, FALSE, DEFAULT_CORETIMES, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/overload", "should edge finding be used to detect an overload?",
         &conshdlrdata->overload, FALSE, DEFAULT_OVERLOAD, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/edgefinding", "should edge-finding be executed?",
         &conshdlrdata->edgefinding, FALSE, DEFAULT_EDGEFINDING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/useadjustedjobs", "should edge-finding be executed?",
         &conshdlrdata->useadjustedjobs, FALSE, DEFAULT_USEADJUSTEDJOBS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/usebinvars", "should the binary representation be used?",
         &conshdlrdata->usebinvars, FALSE, DEFAULT_USEBINVARS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/localcuts", "should cuts be added only locally?",
         &conshdlrdata->localcuts, FALSE, DEFAULT_LOCALCUTS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/usecovercuts", "should covering cuts be added every node?",
         &conshdlrdata->usecovercuts, FALSE, DEFAULT_USECOVERCUTS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/cutsasconss",
         "should the Spread constraint create cuts as knapsack constraints?",
         &conshdlrdata->cutsasconss, FALSE, DEFAULT_CUTSASCONSS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/sepaold",
         "shall old sepa algo be applied?",
         &conshdlrdata->sepaold, FALSE, DEFAULT_SEPAOLD, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/fillbranchcands", "should branching candidates be added to storage?",
         &conshdlrdata->fillbranchcands, FALSE, DEFAULT_FILLBRANCHCANDS, NULL, NULL) );

   /* presolving parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/dualpresolve", "should dual presolving be applied?",
         &conshdlrdata->dualpresolve, FALSE, DEFAULT_DUALPRESOLVE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/coeftightening", "should coefficient tightening be applied?",
         &conshdlrdata->coeftightening, FALSE, DEFAULT_COEFTIGHTENING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/normalize", "should demands and capacity be normalized?",
         &conshdlrdata->normalize, FALSE, DEFAULT_NORMALIZE, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip,
         "constraints/"CONSHDLR_NAME"/maxnodes",
         "number of branch-and-bound nodes to solve an independent Spread constraint (-1: no limit)?",
         &conshdlrdata->maxnodes, FALSE, DEFAULT_MAXNODES, -1LL, SCIP_LONGINT_MAX, NULL, NULL) );

   /* conflict analysis parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/usebdwidening", "should bound widening be used during the conflict analysis?",
         &conshdlrdata->usebdwidening, FALSE, DEFAULT_USEBDWIDENING, NULL, NULL) );

   return SCIP_OKAY;
}

//---------------------------------------------------------------------------------------------------------------------
 
/** creates and captures an Spread constraint in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsSpread(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsSpread() for information about the basic constraint flag configuration*/
SCIP_RETCODE SCIPcreateConsBasicSpread(
   SCIP*                 scip,               //< SCIP data structure 
   SCIP_VAR**            vars,               //< array with variables of constraint entries 
   SCIP_CONS**           cons,               //< pointer to hold the created constraint 
   const char*           name,               //< name of constraint 
   int                   nvars              //< number of variables in the constraint 
   )
{
   SCIP_CALL( SCIPcreateConsSpread(scip, vars, cons, name, nvars,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}


