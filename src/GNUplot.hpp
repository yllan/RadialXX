/*------------------------------------------------------------------------
 *  Copyright (C) 2008  Luis M. de la Cruz
 *
 *  This file is part of TUNA::RBF
 *
 *  TUNA::RBF is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  TUNA::RBF is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ------------------------------------------------------------------------*/


#ifdef WITH_GNUPLOT

#ifndef _GNUPLOT_HPP_
#define _GNUPLOT_HPP_

/*!
 *  class for plotting the results in execution time using GNUplot.
 */
class GNUplot {
public:

/*!
 *  The constructor jsut define the "pipe" to controll GNUplot.
 */
    GNUplot() {
   gnuplotpipe = popen("gnuplot","w");
   if (!gnuplotpipe) {
       throw("Gnuplot not found !");
   }
    }
/*!
 *  Free the "pipe".
 */
    ~GNUplot() {
   fprintf(gnuplotpipe,"exit\n");
   pclose(gnuplotpipe);
    }

/*!
 *  operator overloading to send commands to GNUplot.
 */
    void operator()(const std::string& command) {
   fprintf(gnuplotpipe,"%s\n",command.c_str());
   fflush(gnuplotpipe);   
    }

protected:
    FILE *gnuplotpipe; ///< The "pipe" variable. 
};

#endif // _GNUPLOT_H_

#endif // WITH_GNUPLOT
