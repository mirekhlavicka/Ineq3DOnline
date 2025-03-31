using MeshData;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Web;

namespace Ineq3DOnline.Models
{
    public class IneqMeshTools
    {
        const double minQuality = 0.25d;
        public static void CheckQuality(IneqMesh ineqMesh)
        {
            if (ineqMesh == null)
                return;

            int c = 0;
            Dictionary<int, int> counts = new Dictionary<int, int>();
            do
            {
                c = ineqMesh.CheckQuality(minQuality, false);

                if (!counts.Keys.Contains(c))
                    counts[c] = 1;
                else
                    counts[c]++;


            }
            while (c != 0 && counts[c] < 3);

            if (c != 0)
            {
                c = ineqMesh.CheckQuality(minQuality, true);
            }

            ineqMesh.CheckTopology();
            //ineqMesh.DeleteLonelyPoints();

            ineqMesh.Jiggle(3);

            return;
        }

        public static void CheckBoundaryQuality(IneqMesh ineqMesh)
        {
            if (ineqMesh == null)
                return;

            int c = 0;
            Dictionary<int, int> counts = new Dictionary<int, int>();
            do
            {
                c = ineqMesh.CheckBoundaryQuality(minQuality, false);

                if (!counts.Keys.Contains(c))
                    counts[c] = 1;
                else
                    counts[c]++;

            }
            while (c != 0 && counts[c] < 3);

            //if (c != 0)
            //{
            //    c = ineqMesh.CheckBoundaryQuality(minQuality, true);
            //}

            ineqMesh.CheckTopology();
            //ineqMesh.DeleteLonelyPoints();

            ineqMesh.Jiggle(3, false);

            return;
        }

        public static void CheckCurvatureQuality(IneqMesh ineqMesh)
        {
            ineqMesh.CheckCurvatureQuality();
        }
    }
}