using System;
using MeshData;


namespace Ineq3DOnline
{
    using FuncXYZ = Func<double, double, double, double>;
    public static class IneqLib
    {
        public static IneqTree Ball(double x0, double y0, double z0, double r, double p = 2, bool outside = false)
        {
            if (outside)
            {
                return new IneqTree((x, y, z) =>
                    -Math.Pow(Math.Abs(x - x0), p) -
                    Math.Pow(Math.Abs(y - y0), p) -
                    Math.Pow(Math.Abs(z - z0), p) +
                    Math.Pow(r, p));
            }
            else
            {
                return new IneqTree((x, y, z) =>
                    Math.Pow(
                        Math.Pow(Math.Abs(x - x0), p) +
                        Math.Pow(Math.Abs(y - y0), p) +
                        Math.Pow(Math.Abs(z - z0), p), 1/p) -
                    r);
            }
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

        public static IneqTree Cylinder(double x1, double y1, double z1, double x2, double y2, double z2, double r, double d, bool outside = false)
        {
            double nx, ny, nz, n;

            nx = x2 - x1;
            ny = y2 - y1;
            nz = z2 - z1;

            n = Math.Sqrt(nx * nx + ny * ny + nz * nz);
            nx /= n;
            ny /= n;
            nz /= n;

            if (outside)
            {
                return new IneqTree((x, y, z) =>
                {
                    double vp = (x - x1) * nx + (y - y1) * ny + (z - z1) * nz;
                    return
                        -Math.Pow((x - x1) - vp * nx, 2) -
                        Math.Pow((y - y1) - vp * ny, 2) -
                        Math.Pow((z - z1) - vp * nz, 2) +
                        r * r;
                }) |
                    new IneqTree((x, y, z) =>
                    {
                        double vp = (x - (x1 + x2) / 2) * nx + (y - (y1 + y2) / 2) * ny + (z - (z1 + z2) / 2) * nz;
                        return -Math.Abs(vp) + d;
                    });

            }
            else
            {
                return new IneqTree((x, y, z) =>
                {
                    double vp = (x - x1) * nx + (y - y1) * ny + (z - z1) * nz;
                    return
                        Math.Sqrt(Math.Pow((x - x1) - vp * nx, 2) +
                        Math.Pow((y - y1) - vp * ny, 2) +
                        Math.Pow((z - z1) - vp * nz, 2)) -
                        r;
                }) &
                    new IneqTree((x, y, z) =>
                    {
                        double vp = (x - (x1 + x2) / 2) * nx + (y - (y1 + y2) / 2) * ny + (z - (z1 + z2) / 2) * nz;
                        return Math.Abs(vp) - d;
                    });

            }
        }

        public static IneqTree Balls(int count, double R, double r, double p = 2)
        {
            IneqTree res = new IneqTree();

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
            IneqTree res = new IneqTree();

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
            IneqTree res = new IneqTree();

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
            IneqTree res = new IneqTree();

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
            IneqTree res = new IneqTree();

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
            IneqTree res = new IneqTree();

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
            IneqTree res = new IneqTree();

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
            IneqTree res = new IneqTree();

            for (int i = 0; i < count; i++)
            {
                double x0 = 0.75d * Math.Cos(i * 2 * Math.PI / count);
                double y0 = 0.75d * Math.Sin(i * 2 * Math.PI / count);
                double z0 = 0;

                res = res & ((x, y, z) => (x - x0) * x0 + (y - y0) * y0 + (z - z0) * 0.5);
            }
            return res;
        }

        public static void Gradient(Func<double, double, double, double> f, Point P, out double nx, out double ny, out double nz, double d = 1e-4)
        {
            double w = f(P.X, P.Y, P.Z);

            double wx = f(P.X + d, P.Y, P.Z);
            double wy = f(P.X, P.Y + d, P.Z);
            double wz = f(P.X, P.Y, P.Z + d);
            nx = (wx - w) / d;
            ny = (wy - w) / d;
            nz = (wz - w) / d;

            double n = Math.Sqrt(nx * nx + ny * ny + nz * nz);
            if (n != 0)
            {
                nx = nx / n; ny = ny / n; nz = nz / n;
            }
        }

        public static void ProjectToSurface(this IneqMesh mesh, Point P, int ineqNumber)
        {
            mesh.ProjectToSurface(P, 100, ineqNumber, false);
        }

        public static bool ProjectToSurface(Point P, FuncXYZ f, double D = 0.1d, double precision = 100)
        {
            double n1, n2, n3, n;
            double dx, dy, dz;
            double w, wx, wy, wz;
            dx = D / 10000000d;
            dy = D / 10000000d;
            dz = D / 10000000d;

            w = f(P.X, P.Y, P.Z);

            wx = f(P.X + dx, P.Y, P.Z);
            wy = f(P.X, P.Y + dy, P.Z);
            wz = f(P.X, P.Y, P.Z + dz);
            n1 = (wx - w) / dx;
            n2 = (wy - w) / dy;
            n3 = (wz - w) / dz;

            n = Math.Sqrt(n1 * n1 + n2 * n2 + n3 * n3);
            if (n != 0)
            {
                n1 = n1 / n; n2 = n2 / n; n3 = n3 / n;
            }
            else
            {
                return false;
            }
            Point A, B;
            double w1, tet;

            A = new Point(P.X, P.Y, P.Z);
            tet = -w / n;
            B = new Point(P.X + tet * n1, P.Y + tet * n2, P.Z + tet * n3);

            double stopdist = D / precision;
            double errordist = 2 * D / 3.0d;
            int p = 0;
            int maxItCount = 100;

            if (tet > errordist)
            {
                p = maxItCount;
            }

            while (Math.Abs(tet) > stopdist && p < maxItCount)
            {
                p++;
                w1 = f(B.X, B.Y, B.Z);
                if (w1 == w || tet > errordist)
                {
                    p = maxItCount;
                    break;
                }
                tet = -w1 * tet / (w1 - w);
                A.MoveTo(B, false);
                w = w1;
                B.MoveTo(B.X + tet * n1, B.Y + tet * n2, B.Z + tet * n3, false);
            }
            //if (p < maxItCount && Math.Sqrt((P.X - B.X) * (P.X - B.X) + (P.Y - B.Y) * (P.Y - B.Y) + (P.Z - B.Z) * (P.Z - B.Z)) < errordist)
            //{
                return P.MoveTo(B, false);
            //}
            //else
            //    return false;
        }


        public static double[,] ComputeBasis(double nx, double ny, double nz, char axis = 'x')
        {
            // The input vector (nx, ny, nz) is assumed to be a unit vector.
            double[] n = { nx, ny, nz };

            // Choose an arbitrary vector different from n
            double[] t = Math.Abs(nz) < Math.Abs(nx) ? new double[] { 0, 0, 1 } : new double[] { 1, 0, 0 };

            // Compute first perpendicular vector b1 = t × n
            double[] b1 = new double[]
            {
                t[1] * n[2] - t[2] * n[1],
                t[2] * n[0] - t[0] * n[2],
                t[0] * n[1] - t[1] * n[0]
            };

            // Normalize b1
            double normB1 = Math.Sqrt(b1[0] * b1[0] + b1[1] * b1[1] + b1[2] * b1[2]);
            b1[0] /= normB1;
            b1[1] /= normB1;
            b1[2] /= normB1;

            // Compute second perpendicular vector b2 = n × b1
            double[] b2 = new double[]
            {
                n[1] * b1[2] - n[2] * b1[1],
                n[2] * b1[0] - n[0] * b1[2],
                n[0] * b1[1] - n[1] * b1[0]
            };

            /*// Compute determinant of the basis matrix
            double determinant = n[0] * (b1[1] * b2[2] - b1[2] * b2[1])
                               - n[1] * (b1[0] * b2[2] - b1[2] * b2[0])
                               + n[2] * (b1[0] * b2[1] - b1[1] * b2[0]);

            // Ensure right-handed orientation
            if (determinant < 0)
            {
                b2[0] = -b2[0];
                b2[1] = -b2[1];
                b2[2] = -b2[2];
            }*/

            // Return as a 3x3 matrix with rows as vectors
            if (axis == 'x')
            {
                return new double[,] { { n[0], n[1], n[2] }, { b1[0], b1[1], b1[2] }, { b2[0], b2[1], b2[2] } };
            }
            else if (axis == 'y')
            {
                return new double[,] { { b2[0], b2[1], b2[2] }, { n[0], n[1], n[2] }, { b1[0], b1[1], b1[2] } };
            }
            else// (axis == 'z')
            {
                return new double[,] { { b1[0], b1[1], b1[2] }, { b2[0], b2[1], b2[2] }, { n[0], n[1], n[2] } };
            }
        }

        public static void Transform(ref double x, ref double y, ref double z, Point P, double[,] m)
        {
            x = x - P.X;
            y = y - P.Y;
            z = z - P.Z;

            double nx = m[0, 0] * x + m[0, 1] * y + m[0, 2] * z;
            double ny = m[1, 0] * x + m[1, 1] * y + m[1, 2] * z;
            double nz = m[2, 0] * x + m[2, 1] * y + m[2, 2] * z;

            x = nx;
            y = ny;
            z = nz;
        }
    }
}