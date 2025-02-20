using MeshData;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Web;
using System.Xml.Linq;

namespace Ineq3DOnline
{
	public class SignedDistance
	{
        private TriangleMeshInterop.TriangleMeshWrapper meshWrapper = null;
        public SignedDistance(string filePath, double d = 2, double x0 = 0, double y0 = 0, double z0 = 0)
		{
            List<double[]> v = null;
            List<int[]> t = null;
            filePath = System.IO.Path.Combine(HttpContext.Current.Server.MapPath("~/MeshSamples"), filePath);

            STLLoader.LoadSTL(filePath, out v, out t);

            if (d > 0)
            {
                double minX = double.MaxValue, minY = double.MaxValue, minZ = double.MaxValue;
                double maxX = -double.MaxValue, maxY = -double.MaxValue, maxZ = -double.MaxValue;

                foreach (var p in v)
                {
                    minX = Math.Min(minX, p[0]);
                    minY = Math.Min(minY, p[1]);
                    minZ = Math.Min(minZ, p[2]);
                    maxX = Math.Max(maxX, p[0]);
                    maxY = Math.Max(maxY, p[1]);
                    maxZ = Math.Max(maxZ, p[2]);
                }

                double s = Math.Min(Math.Min(d / (maxX - minX), d / (maxY - minY)), d / (maxZ - minZ));
                double x00 = (maxX + minX) / 2.0d;
                double y00 = (maxY + minY) / 2.0d;
                double z00 = (maxZ + minZ) / 2.0d;

                foreach (var p in v)
                {
                    p[0] = x0 + s * (p[0] - x00);
                    p[1] = y0 + s * (p[1] - y00);
                    p[2] = z0 + s * (p[2] - z00);
                }
            }

            meshWrapper = new TriangleMeshInterop.TriangleMeshWrapper(v, t);
        }

        public double From(double x, double y, double z)
        {
            double[] result = meshWrapper.SignedDistance(x, y, z);
            return result[0];
        }
    }
}