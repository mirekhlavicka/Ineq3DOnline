using System;
using System.Collections.Generic;
using System.Linq;
using System.Web;
using MeshData;

namespace Ineq3DOnline.PlatonicSolids
{
    public struct Vec3
    {
        public double X, Y, Z;
        public Vec3(double x, double y, double z) { X = x; Y = y; Z = z; }

        public static Vec3 operator +(Vec3 a, Vec3 b) => new Vec3(a.X + b.X, a.Y + b.Y, a.Z + b.Z);
        public static Vec3 operator /(Vec3 a, double d) => new Vec3(a.X / d, a.Y / d, a.Z / d);
        public static Vec3 operator *(double f, Vec3 a) => new Vec3(f * a.X, f * a.Y, f * a.Z);

        public double Length() => Math.Sqrt(X * X + Y * Y + Z * Z);
        public Vec3 Normalize(double r = 1.0d)
        {
            var len = Length();
            return new Vec3(r * X / len, r * Y / len, r * Z / len);
        }

        public static double Dot(Vec3 a, Vec3 b) => a.X * b.X + a.Y * b.Y + a.Z * b.Z;
        public static Vec3 Cross(Vec3 a, Vec3 b) =>
            new Vec3(a.Y * b.Z - a.Z * b.Y, a.Z * b.X - a.X * b.Z, a.X * b.Y - a.Y * b.X);

        public override string ToString() => $"{{ {X:G15}, {Y:G15}, {Z:G15} }}";
    }

    public class Polyhedron
    {
        public IReadOnlyList<Vec3> Vertices { get; }
        public IReadOnlyList<int[]> Faces { get; }
        public IReadOnlyList<(int, int)> Edges { get; }

        private Polyhedron(List<Vec3> vertices, List<int[]> faces)
        {
            Vertices = vertices;
            Faces = faces;

            // build unique unordered edges (min,max)
            var edgeSet = new HashSet<(int, int)>();
            foreach (var face in faces)
            {
                for (int i = 0; i < face.Length; i++)
                {
                    int a = face[i];
                    int b = face[(i + 1) % face.Length];
                    if (a > b) (a, b) = (b, a);
                    edgeSet.Add((a, b));
                }
            }
            Edges = edgeSet.ToList();
        }

        public IneqTree ToIneqTree()
        { 
            var res = new IneqTree();

            foreach (var face in Faces)
            { 
                var vx = face.Sum(i => Vertices[i].X) / face.Length;
                var vy = face.Sum(i => Vertices[i].Y) / face.Length;
                var vz = face.Sum(i => Vertices[i].Z) / face.Length;

                res = res & ((x, y ,z) => (x - vx) * vx + (y - vy) * vy + (z - vz) * vz) ;
            }

            return res;
        }

        public static Polyhedron CreateIcosahedron(double r = 1.0d)
        {
            double phi = (1.0 + Math.Sqrt(5.0)) / 2.0;

            var verts = new List<Vec3>()
        {
            new Vec3(-1,  phi,  0), new Vec3( 1,  phi,  0),
            new Vec3(-1, -phi,  0), new Vec3( 1, -phi,  0),
            new Vec3( 0, -1,  phi), new Vec3( 0,  1,  phi),
            new Vec3( 0, -1, -phi), new Vec3( 0,  1, -phi),
            new Vec3( phi,  0, -1), new Vec3( phi,  0,  1),
            new Vec3(-phi,  0, -1), new Vec3(-phi,  0,  1)
        }.Select(v => v.Normalize(r)).ToList();

            var faces = new List<int[]>
        {
            new[]{0,11,5}, new[]{0,5,1}, new[]{0,1,7}, new[]{0,7,10}, new[]{0,10,11},
            new[]{1,5,9}, new[]{5,11,4}, new[]{11,10,2}, new[]{10,7,6}, new[]{7,1,8},
            new[]{3,9,4}, new[]{3,4,2}, new[]{3,2,6}, new[]{3,6,8}, new[]{3,8,9},
            new[]{4,9,5}, new[]{2,4,11}, new[]{6,2,10}, new[]{8,6,7}, new[]{9,8,1}
        };

            return new Polyhedron(verts, faces);
        }

        public static Polyhedron CreateDodecahedron(double r = 1.0d)
        {
            var ico = CreateIcosahedron();

            // Step 1: face centroids of icosahedron => dodecahedron vertices
            var dodeVerts = new List<Vec3>();
            for (int f = 0; f < ico.Faces.Count; f++)
            {
                var face = ico.Faces[f];
                var centroid = (ico.Vertices[face[0]] + ico.Vertices[face[1]] + ico.Vertices[face[2]]) / 3.0;
                dodeVerts.Add(centroid.Normalize(r));
            }

            // Step 2: for each icosahedron vertex collect adjacent faces and order them
            var dodeFaces = new List<int[]>();
            for (int v = 0; v < ico.Vertices.Count; v++)
            {
                // collect incident face indices
                var incident = new List<int>();
                for (int fi = 0; fi < ico.Faces.Count; fi++)
                {
                    var face = ico.Faces[fi];
                    if (face.Contains(v)) incident.Add(fi);
                }

                if (incident.Count != 5)
                    throw new Exception($"Icosahedron vertex {v} belongs to {incident.Count} faces (expected 5).");

                // ordering around vector ico.Vertices[v] using same frame from original generator
                var vert = ico.Vertices[v];
                Vec3 a = new Vec3(1.0, 0.0, 0.0);
                if (Math.Abs(Vec3.Dot(a, vert)) > 0.9) a = new Vec3(0.0, 1.0, 0.0);

                var u = Vec3.Cross(a, vert).Normalize();
                var w = Vec3.Cross(vert, u).Normalize();

                var angleList = incident
                    .Select(fi =>
                    {
                        var fpt = dodeVerts[fi];
                        double x = Vec3.Dot(fpt, u);
                        double y = Vec3.Dot(fpt, w);
                        double ang = Math.Atan2(y, x);
                        return new { fi, ang };
                    })
                    .OrderBy(p => p.ang)
                    .Select(p => p.fi)
                    .ToArray();

                dodeFaces.Add(angleList);
            }

            return new Polyhedron(dodeVerts, dodeFaces);
        }
    }
}