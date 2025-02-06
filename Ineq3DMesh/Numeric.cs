using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MeshData
{
    public class Numeric
    {
        public static bool Newton2d(ref double x0, ref double y0, Func<double, double, double> func1, Func<double, double, double> func2, double D, int maxIterations = 100)
        {
            double x, y, f1, f2, f1x, f1y, f2x, f2y, h1 = D, h2 = D;
            double dx, dy;
            double det;
            dx = D / 10000000d;
            dy = D / 10000000d;

            x = x0; y = y0;
            int iteration = 0;
            while (iteration < maxIterations)
            {
                f1 = func1(x, y);
                f2 = func2(x, y);

                f1x = (func1(x + dx, y) - f1) / dx;
                f1y = (func1(x, y + dy) - f1) / dy;
                f2x = (func2(x + dx, y) - f2) / dx;
                f2y = (func2(x, y + dy) - f2) / dy;

                det = f1x * f2y - f1y * f2x;
                if (det == 0)
                    return false;

                h1 = (-f1 * f2y + f1y * f2) / det;
                h2 = (-f1x * f2 + f1 * f2x) / det;
                x += h1; y += h2;

                if (Math.Abs(h1) + Math.Abs(h2) < D / 200d)
                {
                    break;
                }

                iteration++;
            }
            if (Math.Abs(h1) + Math.Abs(h2) < D / 200d)
            {
                x0 = x;
                y0 = y;
                return true;
            }
            else
            {
                return false;
            }
        }

        public static bool Newton3d(ref double x0, ref double y0, ref double z0, Func<double, double, double, double> func1, Func<double, double, double, double> func2, Func<double, double, double, double> func3, double D, int maxIterations = 100)
        {
            double x, y, z, f1, f2, f3, f1x, f1y, f1z, f2x, f2y, f2z, f3x, f3y, f3z, h1 = D, h2 = D, h3 = D;
            double dx, dy, dz;
            double det;
            dx = D / 10000000d;
            dy = D / 10000000d;
            dz = D / 10000000d;

            x = x0; y = y0; z = z0;
            int iteration = 0;
            while (iteration < maxIterations)
            {
                f1 = func1(x, y, z);
                f2 = func2(x, y, z);
                f3 = func3(x, y, z);

                f1x = (func1(x + dx, y, z) - f1) / dx;
                f1y = (func1(x, y + dy, z) - f1) / dy;
                f1z = (func1(x, y, z + dz) - f1) / dz;

                f2x = (func2(x + dx, y, z) - f2) / dx;
                f2y = (func2(x, y + dy, z) - f2) / dy;
                f2z = (func2(x, y, z + dz) - f2) / dz;

                f3x = (func3(x + dx, y, z) - f3) / dx;
                f3y = (func3(x, y + dy, z) - f3) / dy;
                f3z = (func3(x, y, z + dz) - f3) / dz;


                det = Determinant(f1x, f1y, f1z, f2x, f2y, f2z, f3x, f3y, f3z);
                if (det == 0)
                    return false;

                h1 = Determinant(f1, f1y, f1z, f2, f2y, f2z, f3, f3y, f3z) / det;
                h2 = Determinant(f1x, f1, f1z, f2x, f2, f2z, f3x, f3, f3z) / det;
                h3 = Determinant(f1x, f1y, f1, f2x, f2y, f2, f3x, f3y, f3) / det;

                x -= h1; y -= h2; z -= h3;
                if(Math.Abs(h1) + Math.Abs(h2) + Math.Abs(h3) < D / 200d)
                {
                    break;
                }

                iteration++;

            }
            if (Math.Abs(h1) + Math.Abs(h2) + Math.Abs(h3) < D / 200d)
            {
                x0 = x;
                y0 = y;
                z0 = z;
                return true;
            }
            else
            {
                return false;
            }
        }

        public static double Determinant(double a11, double a12, double a13, double a21, double a22, double a23, double a31, double a32, double a33)
        {
            return a11 * a22 * a33 + a21 * a32 * a13 + a12 * a23 * a31 - a13 * a22 * a31 - a12 * a21 * a33 - a23 * a32 * a11;
        }

        public static double Interpolate(
            double f000, double f010, double f001, double f011,
            double f100, double f110, double f101, double f111,
            double x, double y, double z)
        {
            return
                f000 * (1 - x) * (1 - y) * (1 - z) +
                f010 * (1 - x) * y * (1 - z) +
                f001 * (1 - x) * (1 - y) * z +
                f011 * (1 - x) * y * z +

                f100 * x * (1 - y) * (1 - z) +
                f110 * x * y * (1 - z) +
                f101 * x * (1 - y) * z +
                f111 * x * y * z;
        }
    }
}
