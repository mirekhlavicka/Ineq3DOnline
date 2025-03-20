using System;
using System.Collections.Generic;
using System.Linq;
using MeshData;
using Ineq3DOnline;
using static Ineq3DOnline.MyMath;
using static Ineq3DOnline.IneqLib;
using Jace.Operations;

namespace TempNamespace
{
    public class TempClass
    {
        private SignedDistance sdf = new SignedDistance("Henne.stl");

        public IneqTree BackEgg()
        {
            IneqTree res = new IneqTree();

            Point p = new Point(-0.35d, 1.0d, -0.4d);

            sdf.Project(p);
            Gradient(sdf.From, p, out double nx, out double ny, out double nz);

            p.X += 0.05 * nx;
            p.Y += 0.05 * ny;
            p.Z += 0.05 * nz;


            var m = ComputeBasis(nx, ny, nz, 'z');

            res = res | new IneqTree((x, y, z) =>
            {
                Transform(ref x, ref y, ref z, p, m);

                x = x * 10;
                y = y * 10;
                z = z * 10;

                return sqrt(pow(pow(x, 2) + pow(y, 2) + pow(z, 2) / 2.25d, 3)) - pow(x, 2) - pow(y, 2) - pow(z, 2) / 2.25d;
            });
            return res;
        }

        public IneqTree Eggs(int count)
        {
            IneqTree res = new IneqTree();
            Random random = new Random();
            List<Point> points = new List<Point>();

            for (int i = 0; i < count; i++)
            {
                Point p = new Point(
                        random.NextDouble() * 1.6d - 0.8d,
                        random.NextDouble() * 1.1d - 0.8d,
                        random.NextDouble() * 0.1d - 0.7d);

                if (abs(sdf.From(p.X, p.Y, p.Z)) < 0.1d)
                {
                    sdf.Project(p);
                    if (!points.Any(p1 => p1.Distance(p) < 2.3d * 0.1d))
                    {
                        Gradient(sdf.From, p, out double nx, out double ny, out double nz);

                        p.X += 0.05 * nx;
                        p.Y += 0.05 * ny;
                        p.Z += 0.05 * nz;

                        points.Add(p);

                        var m = ComputeBasis(nx, ny, nz, 'x');

                        res = res | new IneqTree((x, y, z) =>
                        {
                            Transform(ref x, ref y, ref z, p, m);

                            x = x * 10;
                            y = y * 10;
                            z = z * 10;

                            return sqrt(pow(pow(x, 2) + pow(y, 2) + pow(z, 2) / 2.25d, 3)) - pow(x, 2) - pow(y, 2) - pow(z, 2) / 2.25d;
                        });
                    }
                }
            }
            return res;
        }

        public IneqMesh GetIneqMesh()
        {
            Func<double, double, double, double> f = sdf.From;
            var b = new IneqTree(f);

            var res = new IneqMesh
            {
                X0 = -1,
                Y0 = -1,
                Z0 = -1,
                X1 = 1,
                Y1 = 1,
                Z1 = 1,
                D = 0.05,
                Boxed = true,
                IneqTree =
                    (b | BackEgg() | Eggs(600)) & !new IneqTree((x, y, z) => z + 0.95)
            };
            res.ProjectToSurfaceSpec[f] = sdf.Project;
            res.PrepareBackgroundMesh = () => res.RefineTetrahedralMeshByTetrahedrons(0, 10, 0.1);

            return res;
        }
    }
}