using System;
using System.Collections.Generic;
using System.Linq;
using System.Web;

namespace Ineq3DOnline
{
    public static class HelixDistance
    {
        public static double Squared(double x, double y, double z, double r = 1, double c = 1, double tolerance = 1e-6, int maxIterations = 100)
        {
            double t = z / c;
            double t0 = z / c -  1 * Math.PI; 

            while (t0 < z / c + 1 * Math.PI)
            {
                if(SquaredDistance(x, y, z, r, c, t0) < SquaredDistance(x, y, z, r, c, t)) 
                { 
                    t = t0;
                }
                t0 += 0.1;
            }

            int iteration = 0;


            while (iteration < maxIterations)
            {
                double fd = FirstDerivative(x, y, z, r, c, t);
                double sd = SecondDerivative(x, y, z, r, c, t);

                if (Math.Abs(sd) < 1e-8)
                {
                    break;
                    //throw new Exception("Second derivative too small, Newton's method failed.");
                }

                double tNext = t - fd / sd;

                if (Math.Abs(tNext - t) < tolerance)
                {
                    t = tNext;
                    break;
                }

                t = tNext;
                iteration++;
            }

            /*if (iteration == maxIterations)
            {
                throw new Exception("Newton's method did not converge within the maximum number of iterations.");
            }*/

            return SquaredDistance(x, y, z, r, c, t);
        }

        private static double SquaredDistance(double x, double y, double z, double r, double c, double t)
        {
            double helixX = r * Math.Cos(t);
            double helixY = r * Math.Sin(t);
            double helixZ = c * t;

            double dx = x - helixX;
            double dy = y - helixY;
            double dz = z - helixZ;

            return dx * dx + dy * dy + dz * dz;
        }

        private static double FirstDerivative(double x, double y, double z, double r, double c, double t)
        {
            return r * x * Math.Sin(t) - r * y * Math.Cos(t) - z + c * t;
        }

        private static double SecondDerivative(double x, double y, double z, double r, double c, double t)
        {
            return r * x * Math.Cos(t) + r * y * Math.Sin(t) + c;
        }
    }
}