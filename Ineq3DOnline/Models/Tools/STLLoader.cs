using System;
using System.Collections.Generic;
using System.IO;
using System.Globalization;

class STLLoader
{
    public static void LoadSTL(string filePath, out List<double[]> vertices, out List<int[]> triangles)
    {
        //vertices.Clear();
        //triangles.Clear();

        vertices = new List<double[]>();
        triangles = new List<int[]>();

        using (var reader = new BinaryReader(File.Open(filePath, FileMode.Open, FileAccess.Read)))
        {
            if (IsBinarySTL(reader))
            {
                ReadBinarySTL(reader, vertices, triangles);
            }
            else
            {
                //reader.BaseStream.Seek(0, SeekOrigin.Begin); // Reset stream for text reading
                reader.Close();
                ReadAsciiSTL(filePath, vertices, triangles);
            }
        }
    }

    private static bool IsBinarySTL(BinaryReader reader)
    {
        reader.BaseStream.Seek(80, SeekOrigin.Begin); // Skip 80-byte header
        int numTriangles = reader.ReadInt32();
        long expectedSize = 80 + 4 + (numTriangles * (50));
        return reader.BaseStream.Length == expectedSize;
    }

    private static void ReadBinarySTL(BinaryReader reader, List<double[]> vertices, List<int[]> triangles)
    {
        reader.BaseStream.Seek(80, SeekOrigin.Begin);
        int numTriangles = reader.ReadInt32();
        Dictionary<(float, float, float), int> vertexMap = new Dictionary<(float, float, float), int>();

        for (int i = 0; i < numTriangles; i++)
        {
            reader.ReadBytes(12); // Skip normal vector
            int[] indices = new int[3];

            for (int j = 0; j < 3; j++)
            {
                float x = reader.ReadSingle();
                float y = reader.ReadSingle();
                float z = reader.ReadSingle();
                var key = (x, y, z);

                if (!vertexMap.TryGetValue(key, out int index))
                {
                    index = vertices.Count;
                    vertices.Add(new double[] { x, y, z });
                    vertexMap[key] = index;
                }

                indices[j] = index;
            }

            triangles.Add(indices);
            reader.ReadBytes(2); // Skip attribute byte count
        }
    }

    private static void ReadAsciiSTL(string filePath, List<double[]> vertices, List<int[]> triangles)
    {
        Dictionary<(double, double, double), int> vertexMap = new Dictionary<(double, double, double), int>();
        using (var reader = new StreamReader(filePath))
        {
            string line;
            int[] indices = new int[3];
            int vertexIndex = 0;
            char[] separators = new char[] { ' '};

            while ((line = reader.ReadLine()) != null)
            {
                line = line.Trim();
                if (line.StartsWith("vertex"))
                {                                   
                    string[] parts = line.Split(separators, StringSplitOptions.RemoveEmptyEntries);
                    double x = double.Parse(parts[1], CultureInfo.InvariantCulture);
                    double y = double.Parse(parts[2], CultureInfo.InvariantCulture);
                    double z = double.Parse(parts[3], CultureInfo.InvariantCulture);
                    var key = (x, y, z);

                    if (!vertexMap.TryGetValue(key, out int index))
                    {
                        index = vertices.Count;
                        vertices.Add(new double[] { x, y, z });
                        vertexMap[key] = index;
                    }

                    indices[vertexIndex++] = index;
                    if (vertexIndex == 3)
                    {
                        triangles.Add(indices);
                        indices = new int[3];
                        vertexIndex = 0;
                    }
                }
            }
        }
    }
}

