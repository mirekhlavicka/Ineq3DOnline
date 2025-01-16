using System;
using System.Collections.Generic;
using System.Linq;
using System.Web;

namespace Ineq3DOnline
{
    public static class HelixDistance
    {
        public static double Squared(
             double x, double y, double z, double r = 1, double c = 1, double tolerance = 1e-6, int maxIterations = 100)
        {
            double t = z / c; // Initial guess
            int iteration = 0;

            while (iteration < maxIterations)
            {
                // Compute first and second derivatives of the squared distance function
                double fPrime = ComputeFirstDerivative(x, y, z, r, c, t);
                double fDoublePrime = ComputeSecondDerivative(x, y, z, r, c, t);

                // Avoid division by zero
                if (Math.Abs(fDoublePrime) < 1e-8)
                {
                    break;
                    //throw new Exception("Second derivative too small, Newton's method failed.");
                }

                // Newton's method update
                double tNext = t - fPrime / fDoublePrime;

                // Check for convergence
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

            // Compute the minimum squared distance
            //var res = ComputeSquaredDistance(x, y, z, r, c, t); 
            return ComputeSquaredDistance(x, y, z, r, c, t);
        }

        private static double ComputeSquaredDistance(double x, double y, double z, double r, double c, double t)
        {
            double helixX = r * Math.Cos(t);
            double helixY = r * Math.Sin(t);
            double helixZ = c * t;

            double dx = x - helixX;
            double dy = y - helixY;
            double dz = z - helixZ;

            return dx * dx + dy * dy + dz * dz;
        }

        private static double ComputeFirstDerivative(double x, double y, double z, double r, double c, double t)
        {
            return r * x * Math.Sin(t) - r * y * Math.Cos(t) - z + c * t;
        }

        private static double ComputeSecondDerivative(double x, double y, double z, double r, double c, double t)
        {
            return r * x * Math.Cos(t) + r * y * Math.Sin(t) + c;
        }
    }
}