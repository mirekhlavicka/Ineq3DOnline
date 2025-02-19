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
        public SignedDistance(string filePath)
		{
            List<double[]> v = null;
            List<int[]> t = null;
            filePath = System.IO.Path.Combine(HttpContext.Current.Server.MapPath("~/MeshSamples"), filePath);

            STLLoader.LoadSTL(filePath, out v, out t);

            meshWrapper = new TriangleMeshInterop.TriangleMeshWrapper(v, t);
        }

        public double From(double x, double y, double z)
        {
            double[] result = meshWrapper.SignedDistance(x, y, z);
            return result[0];
        }
    }
}