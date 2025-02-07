using System;
using MeshData;


namespace Ineq3DOnline
{
    public class IneqLib
    {
        public static IneqTree Ball(double x0, double y0, double z0, double r, double p = 2)
        {
            return new IneqTree((x, y, z) => 
                Math.Pow(Math.Abs(x - x0), p) +
                Math.Pow(Math.Abs(y - y0), p) +
                Math.Pow(Math.Abs(z - z0), p) -
                Math.Pow(r, p));
        }

        public static IneqTree Torus(double x0, double y0, double z0, double r, double R, double nx, double ny, double nz, bool outside = false)
        {

            return new IneqTree((x, y, z) =>
            {
                x = x - x0;
                y = y - y0;
                z = z - z0;
                double vp = x * nx + y * ny + z * nz;
                double v1 = x - vp * nx;
                double v2 = y - vp * ny;
                double v3 = z - vp * nz;
                double v = Math.Sqrt(v1 * v1 + v2 * v2 + v3 * v3);

                if (outside)
                {
                    return -Math.Pow(R - v, 2) - vp * vp + r * r;
                }
                else
                {
                    return Math.Pow(R - v, 2) + vp * vp - r * r;
                }
            });
        }

        public static IneqTree Cylinder(double x1, double y1, double z1, double x2, double y2, double z2, double r, double d)
        {
            double nx, ny, nz, n;

            nx = x2 - x1;
            ny = y2 - y1;
            nz = z2 - z1;

            n = Math.Sqrt(nx * nx + ny * ny + nz * nz);
            nx /= n;
            ny /= n;
            nz /= n;

            return new IneqTree((x, y, z) =>
            {
                double vp = (x - x1) * nx + (y - y1) * ny + (z - z1) * nz;
                return
                    Math.Pow((x - x1) - vp * nx, 2) +
                    Math.Pow((y - y1) - vp * ny, 2) +
                    Math.Pow((z - z1) - vp * nz, 2) -
                    r * r;
            }) &
                new IneqTree((x, y, z) =>
                {
                    double vp = (x - (x1 + x2) / 2) * nx + (y - (y1 + y2) / 2) * ny + (z - (z1 + z2) / 2) * nz;
                    return Math.Abs(vp) - d;
                });
        }

        public static IneqTree Balls(int count, double R, double r, double p = 2)
        {
            IneqTree res = new IneqTree((x, y, z) => 1);

            for (int i = 0; i < count; i++)
            {
                double x0 = R * Math.Cos(i * 2 * Math.PI / count);
                double y0 = R * Math.Sin(i * 2 * Math.PI / count);

                res = res | ((x, y, z) => Math.Pow(Math.Abs(x - x0), p) +
                    Math.Pow(Math.Abs(y - y0), p) +
                    Math.Pow(Math.Abs(z), p) -
                    Math.Pow(r, p));
            }

            return res;
        }

        public static IneqTree Toruses(int count, double z0, double z1, double r0, double r1, double R0, double R1)
        {
            IneqTree res = new IneqTree((x, y, z) => 1);

            for (int i = 0; i < count; i++)
            {
                double z = z0 + i * (z1 - z0) / (count - 1);
                double r = r0 + i * (r1 - r0) / (count - 1);
                double R = R0 + i * (R1 - R0) / (count - 1);

                res = res | Torus(0, 0, z, r, R, 0, 0, 1);
            }

            return res;
        }

        public static IneqTree Toruses1(int count, double RR, double R, double r)
        {
            IneqTree res = new IneqTree((x, y, z) => 1);

            for (int i = 0; i < count; i++)
            {
                double rx = Math.Cos(i * 2 * Math.PI / count);
                double ry = Math.Sin(i * 2 * Math.PI / count);


                double x0 = RR * rx;
                double y0 = RR * ry;

                res = res | Torus(x0, y0, 0, r, R, -ry, rx, 0);
            }

            return res;
        }

        public static IneqTree AntiToruses(int count, double RR, double R, double r)
        {
            IneqTree res = new IneqTree((x, y, z) => -1);

            for (int i = 0; i < count; i++)
            {
                double rx = Math.Cos(i * 2 * Math.PI / count);
                double ry = Math.Sin(i * 2 * Math.PI / count);


                double x0 = RR * rx;
                double y0 = RR * ry;

                res = res & Torus(x0, y0, 0, r, R, -ry, rx, 0, true);
            }

            return res;
        }

        public static IneqTree Cylinders(int count, double R, double r)
        {
            IneqTree res = new IneqTree((x, y, z) => 1);

            for (int i = 0; i < count; i++)
            {
                double x0 = R * Math.Cos(i * 2 * Math.PI / count);
                double y0 = R * Math.Sin(i * 2 * Math.PI / count);

                res = res | ((x, y, z) => (x - x0) * (x - x0) + (y - y0) * (y - y0) - r * r);
            }

            return res;
        }

        public static IneqTree Cylinders1(int count, double R, double r)
        {
            IneqTree res = new IneqTree((x, y, z) => 1);

            for (int i = 0; i < count; i++)
            {
                double x0 = Math.Cos(i * 2 * Math.PI / count);
                double y0 = Math.Sin(i * 2 * Math.PI / count);

                res = res | (IneqTree)(((x, y, z) =>
                {
                    double vp = x * x0 + y * y0;
                    return vp < 0 ? 1 : (x - vp * x0) * (x - vp * x0) + (y - vp * y0) * (y - vp * y0) + z * z - r * r;
                }) &
                (IneqTree)
                    ((x, y, z) =>
                    {
                        double vp = x * x0 + y * y0;
                        return -vp + R;
                    }
                ));
            }

            return res;
        }

        public static IneqTree Planes(int count)
        {
            IneqTree res = new IneqTree((x, y, z) => -1);

            for (int i = 0; i < count; i++)
            {
                double x0 = 0.75d * Math.Cos(i * 2 * Math.PI / count);
                double y0 = 0.75d * Math.Sin(i * 2 * Math.PI / count);
                double z0 = 0;

                res = res & ((x, y, z) => (x - x0) * x0 + (y - y0) * y0 + (z - z0) * z0);
            }

            return res;
        }

        public static IneqTree Planes1(int count)
        {
            IneqTree res = new IneqTree((x, y, z) => -1);

            for (int i = 0; i < count; i++)
            {
                double x0 = 0.75d * Math.Cos(i * 2 * Math.PI / count);
                double y0 = 0.75d * Math.Sin(i * 2 * Math.PI / count);
                double z0 = 0;

                res = res & ((x, y, z) => (x - x0) * x0 + (y - y0) * y0 + (z - z0) * 0.5);
            }
            return res;
        }

    }
}