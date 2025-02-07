using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MeshData;
using Jace;
using System.Web;
using System.Xml.Linq;
using System.Web.Helpers;

namespace Ineq3DOnline
{
    public class MeshSamples
    {
        private Dictionary<string, object> samples = new Dictionary<string, object>();

        public IEnumerable<string> Samples
        {
            get { return samples.Keys; }
        }

        public MeshSamples()
        {
            samples["Half ball and balls"] = 
            new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -1.0,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 1.0,
                D = 0.1d,
                Boxed = false,
                IneqTree =
                    (
                        ((IneqTree)((x, y, z) => x * x + y * y + z * z - 1)) &
                        (
                            (IneqTree)((x, y, z) => -x ) |
                            ((x, y, z) => -y
                        )) & 
                        ((x, y, z) => -z
                        ) | 
                        ((x, y, z) => x*x + y*y + z*z - 0.25) |
                        ((x, y, z) => x*x + y*y + (z+0.5)*(z+0.5) - 0.15)                        
                    )
            };

            samples["Cylinders, ball, planes"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -1.0,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 1.0,
                D = 0.1d,
                Boxed = false,
                IneqTree =
                    (
                        ((IneqTree)
                            ((IneqTree)((x, y, z) => x * x + y * y + z * z - 1) &
                            ((x, y, z) => x - 0.25) &
                            ((x, y, z) => y - 0.25)) |
                            ((x, y, z) => y * y + z * z - 0.2) |
                            ((x, y, z) => x * x + z * z - 0.2)
                        ) &
                            ((x, y, z) => Math.Abs(x) - 1.0) &
                            ((x, y, z) => Math.Abs(y) - 1.0)

                    )
            };

            samples["4-Ball, cylinder, ball, planes"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -1.0,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 1.0,
                D = 0.1d,
                Boxed = false,
                IneqTree =
                    (
                        ((IneqTree)
                            ((x, y, z) => x * x * x * x + y * y * y * y + z * z * z * z - 1) &
                            ((x, y, z) => x - y - z - 0.2d)
                        ) |
                        ((x, y, z) => x * x + y * y + z * z - 0.2)
                    ) &
                    ((x, y, z) => -x * x - y * y + 0.1) &
                    ((x, y, z) => -z - 0.5)            
            };

            samples["Cylinders"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -1.0,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 1.0,
                D = 0.1d,
                Boxed = true,
                IneqTree =
                    (
                        (IneqTree)
                            ((x, y, z) => y * y + z * z - 0.15) |
                            ((x, y, z) => x * x + z * z - 0.15) |
                            ((x, y, z) => x * x + y * y - 0.15)
                    ) &
                    (
                        (IneqTree)
                        ((x, y, z) => -x * x - y * y + 0.05) &
                        ((x, y, z) => -x * x - z * z + 0.05) &
                        ((x, y, z) => -z * z - y * y + 0.05)
                    ) &
                    ((x, y, z) => x)
            };

            samples["Sin 1"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -1.0,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 1.0,
                D = 0.1d,
                Boxed = false,
                IneqTree =
                    ((IneqTree)((x, y, z) => Math.Cos(5 * x) + Math.Cos(5 * y) + Math.Cos(5 * z) - 0.5) &
                    ((x, y, z) => x * x + y * y + z * z - 1) &
                    ((x, y, z) => x - 0.1))
            };

            samples["Sin 2"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -1.0,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 1.0,
                D = 0.05d,
                Boxed = false,
                IneqTree =
                    (
                        (IneqTree)
                            ((x, y, z) => x * x + y * y + z * z - ((Math.Cos(9 * x) + Math.Cos(9 * y) + Math.Cos(9 * z)) / 6.0 + 0.6)) |
                            ((x, y, z) => x * x + y * y + z * z - ((Math.Sin(9 * x) + Math.Sin(9 * y) + Math.Sin(9 * z)) / 6.0 + 0.6))
                    ) &
                    (IneqTree)((x, y, z) => x - 0.1)
            };

            samples["Balls"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -0.25,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 0.25,
                D = 0.075d,
                Boxed = false,
                IneqTree = (IneqTreeBalls(8, 0.5d, 0.267d))
                & ((x, y, z) => z - 0.1)
            };

            samples["Planes"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -0.25,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 0.25,
                D = 0.075d,
                Boxed = false,
                IneqTree = IneqTreePlanes(12) &
                ((x, y, z) => Math.Abs(z) - 0.25) &
                ((x, y, z) => -x * x  - y * y  + 0.1) 
            };

            samples["Anuloid, balls"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -0.35,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 0.35,
                D = 0.075d,
                Boxed = false,
                IneqTree = 
                (IneqTree)
                ((x, y, z) => Math.Pow(Math.Sqrt(x * x + y * y) - 0.70, 2)  + z * z - 0.0315) |
                (IneqTreeBalls(6, 0.7d, Math.Sqrt(0.1d)))
            };

            samples["Many cylinders"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -0.25,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 0.25,
                D = 0.075d,
                Boxed = false,
                IneqTree = (IneqTreeCylinders(8, 0.5d,0.253d))
                & ((x, y, z) => Math.Abs(z) - 0.25)
            };

            samples["Cone"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = 0,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 2.5,
                D = 0.1d,
                Boxed = false,
                IneqTree =
                    (
                        (IneqTree)
                            ((x, y, z) => x * x + y * y - z * z / 4) 
                    ) &
                    (IneqTree)((x, y, z) => -z + 0.4)
                     &
                    (IneqTree)((x, y, z) => x * x + y * y + (z * z) / 2 - 2.25 )
            };
                        
            samples["Cone 1"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -1.0,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 1.0,
                D = 0.055d,
                Boxed = true,
                IneqTree =
                    (
                        (IneqTree)
                            ((x, y, z) => (x - 0.12) * (x - 0.12) + (y + 0.171) * (y + 0.171) - (z - 0.037) * (z - 0.037) / 4)
                    ) &
                    (IneqTree)((x, y, z) => z - 0.037)
                     &
                    (IneqTree)((x, y, z) => (x - 0.12) * (x - 0.12) + (y + 0.171) * (y + 0.171) - 0.15f)
            };

            samples["Ball minus cylinders"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -1.0,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 1.0,
                D = 0.1d,
                Boxed = false,
                IneqTree =
                    ((IneqTree)
                    ((x, y, z) => x*x+y*y+z*z-1) &
                    ((x, y, z) => -x*x-y*y-z*z+25.0/64.0)) &
                    ((x, y, z) => -(x+16.0/25.0)*(x+16.0/25.0)-(y-16.0/25.0)*(y-16.0/25.0)+0.25) &
                    ((x, y, z) => -(y+16.0/25.0)*(y+16.0/25.0)-(z+16.0/25.0)*(z+16.0/25.0)+0.25) &
                    ((x, y, z) => -(z - 16.0 / 25.0) * (z - 16.0 / 25.0) - (x - 16.0 / 25.0) * (x - 16.0 / 25.0) + 0.25)
            };

            samples["Ball in ball"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -1.0,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 1.0,
                D = 0.1d,
                Boxed = false,
                IneqTree =
                    (((IneqTree)
                    ((x, y, z) => x * x + y * y + z * z - 1) &
                    ((x, y, z) => -x * x - y * y + 0.15) &
                    ((x, y, z) => -(x-1) * (x-1) - z * z + 0.25))) |
                    ((IneqTree)((x, y, z) => x * x + y * y + z * z - 0.7))
            };

            samples["Anuloids"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -1,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 1.0,
                D = 0.075d,
                Boxed = false,
                IneqTree =
                (IneqTree)
                ((x, y, z) => (x * x + y * y) * Math.Pow(1 - 0.70 / Math.Sqrt(x * x + y * y), 2) + z * z - 0.0315) |
                ((x, y, z) => (x * x + z * z) * Math.Pow(1 - 0.70 / Math.Sqrt(x * x + z * z), 2) + y * y - 0.0315) |
                ((x, y, z) => (z * z + y * y) * Math.Pow(1 - 0.70 / Math.Sqrt(z * z + y * y), 2) + x * x - 0.0315)
            };

            samples["Anuloids1"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -1.0,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 1.0,
                D = 0.05d,
                Boxed = false,
                IneqTree = IneqTreeAnuloids(6, -0.7d, 0.7d, 0.2d, 0.2d, 0.7d, 0.3d)
                & ((x, y, z) => Math.Abs(z) - 0.8)
            };

            samples["Anuloids2"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -1.0,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 1.0,
                D = 0.045d,
                Boxed = false,
                IneqTree = IneqTreeAnuloids1(8, 0.7d, 0.2d, 0.1d)
                | IneqTreeAnuloid(0, 0, 0, 0.15d, 0.7d, 0, 0, 1)
            };

            samples["Anuloids3"] = new IneqMesh
            {
                X0 = -1.4,
                Y0 = -1.4,
                Z0 = -1.0,
                X1 = 1.4,
                Y1 = 1.4,
                Z1 = 1.0,
                D = 0.065d,
                Boxed = false,
                IneqTree = IneqTreeAnuloids1(16, 0.7d, 0.4d, 0.2d)
                //| IneqTreeAnuloid(0, 0, 0, 0.15d, 0.7d, 0, 0, 1)
            };

            samples["AntiAnuloids"] = new IneqMesh
            {
                X0 = -1.4,
                Y0 = -1.4,
                Z0 = -1.0,
                X1 = 1.4,
                Y1 = 1.4,
                Z1 = 1.0,
                D = 0.065d,
                Boxed = false,
                IneqTree = IneqTreeAntiAnuloids(16, 0.7d, 0.4d, 0.1d)
                & IneqTreeAnuloid(0, 0, 0, 0.4d, 0.7d, 0, 0, 1, false)
            };


            samples["Many cylinders 1"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -0.6,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 0.6,
                D = 0.075d,
                Boxed = false,
                IneqTree = ((IneqTreeCylinders1(6, 0.2d, 0.2d)
                & ((x, y, z) => x * x + y * y + z * z - 1))/*
                | ((x, y, z) => x * x + y * y + 800 * z * z * z * z * z * z * z * z - 0.25d)*/
                )
            };

            samples["Viviani"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -1.0,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 1.0,
                D = 0.1d,
                Boxed = false,
                IneqTree = (IneqTree)((x, y, z) => x * x + y * y + z * z - 1) & ((x, y, z) => (x - 0.5) * (x - 0.5) + y * y - 0.25)                
            };

            samples["Sharp cone"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -1,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 1.2,
                D = 0.075d,
                Boxed = false,
                IneqTree = IneqTreePlanes1(6) &
                    //((x, y, z) => z - 1) &
                ((x, y, z) => -z)
                /*&
                ((x, y, z) => -x * x  - y * y  + 0.1)*/
            };
            samples["Balls in 4-ball"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -1.0,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 1.0,
                D = 0.1d,
                Boxed = false,
                IneqTree =
                (
                    (
                        (IneqTree)((x, y, z) => x * x * x * x + y * y * y * y + z * z * z * z - 1) &
                        ((x, y, z) => -x * x - y * y + 0.37) &
                        ((x, y, z) => x + y + z - 0.4)
                    )
                    | ((x, y, z) => (x - 0.35) * (x - 0.35) + y * y + z * z - 0.30)
                    | ((x, y, z) => (x + 0.35) * (x + 0.35) + y * y + z * z - 0.30)
                )
                &
                    ((x, y, z) => y - 0.7)

            };

            samples["Ball in box"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -1.0,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 1.0,
                D = 0.1d,
                Boxed = false,
                IneqTree =
                (
                    (
                        (IneqTree)((x, y, z) => Math.Abs(x) - 0.6) &
                        ((x, y, z) => Math.Abs(y) - 0.6) &
                        ((x, y, z) => Math.Abs(z) - 0.6)
                    )
                    | ((x, y, z) => x * x + y * y + z * z - 0.5)
                )

            };

            samples["SnowWoman"] = new IneqMesh
            {
                X0 = -0.6,
                Y0 = -1.0,
                Z0 = -0.6,
                X1 = 0.6,
                Y1 = 1.0,
                Z1 = 0.6,
                D = 0.06,
                Boxed = false,
                IneqTree =
                    (
                        ((IneqTree)((x, y, z) => (y - 0.5) * (y - 0.5) + x * x + z * z - 0.29) |
                        ((x, y, z) => (y - 0) * (y - 0) + x * x + z * z - 0.2) |
                        ((x, y, z) => (y - 0) * (y - 0) + (x + 0.4) * (x + 0.4) + (z - 0.2) * (z - 0.2) - 0.05) |
                        ((x, y, z) => (y - 0) * (y - 0) + (x + 0.4) * (x + 0.4) + (z + 0.2) * (z + 0.2) - 0.05) |
                        ((x, y, z) => (y + 0.5) * (y + 0.5) + x * x + z * z - 0.15) |
                        ((x, y, z) => 2 * (y + 0.5) * (y + 0.5) + (x + 0.4) * (x + 0.4) + 2 * z * z - 0.01) |
                        ((x, y, z) => (x * x + z * z - 0.08))) &
                        ((x, y, z) => (Math.Abs(y) - 1.0))
                    )
            };

            samples["Cuted cylinders"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -1.0,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 1.0,
                D = 0.05d,
                Boxed = true,
                IneqTree =
                    (
                        (
                            (IneqTree)((x, y, z) => x * x + y * y - 0.25) |
                            ((x, y, z) => x * x + z * z - 0.25) |
                            ((x, y, z) => z * z + y * y - 0.25)
                        ) &
                        (IneqTree)((x, y, z) => Math.Pow(Math.Abs(x), 0.5) + Math.Pow(Math.Abs(y), 0.5) + Math.Pow(Math.Abs(z), 0.5) - 1.89)
                    )
            };

            samples["Ball"] = new IneqMesh
            {
                X0 = -1.0,
                Y0 = -1.0,
                Z0 = -1.0,
                X1 = 1.0,
                Y1 = 1.0,
                Z1 = 1.0,
                D = 0.15d,
                Boxed = false,
                IneqTree =
                        ((IneqTree)((x, y, z) => x * x + y * y + z * z - 1))
            };

            samples["Tetra cylinder"] = new IneqMesh
            {
                X0 = -0.75,
                Y0 = -0.75,
                Z0 = -0.75,
                X1 = 2.25,
                Y1 = 2.25,
                Z1 = 2.25,
                D = 0.1d,
                Boxed = false,
                IneqTree =
                        IneqTreeBall(0, 0, 0, 0.55) |
                        IneqTreeBall(1.5, 0, 0, 0.55) |
                        IneqTreeBall(0.75, 1.5 * Math.Sin(2 * Math.PI / 3), 0, 0.55) |
                        IneqTreeBall(2.25 / 3, 1.5 * Math.Sin(2 * Math.PI / 3) / 3, Math.Sqrt(2.25 - Math.Pow(2 * 1.5 * Math.Sin(2 * Math.PI / 3) / 3, 2)), 0.55) |

                        IneqTreeCylinder(0, 0, 0, 1.5, 0, 0, 0.2, 1.5) |
                        IneqTreeCylinder(0, 0, 0, 0.75, 1.5 * Math.Sin(2 * Math.PI / 3), 0, 0.2, 1.5) |
                        IneqTreeCylinder(0, 0, 0, 2.25 / 3, 1.5 * Math.Sin(2 * Math.PI / 3) / 3, Math.Sqrt(2.25 - Math.Pow(2 * 1.5 * Math.Sin(2 * Math.PI / 3) / 3, 2)), 0.2, 1.5) |

                        IneqTreeCylinder(1.5, 0, 0, 0.75, 1.5 * Math.Sin(2 * Math.PI / 3), 0, 0.2, 1.5) |
                        IneqTreeCylinder(1.5, 0, 0, 2.25 / 3, 1.5 * Math.Sin(2 * Math.PI / 3) / 3, Math.Sqrt(2.25 - Math.Pow(2 * 1.5 * Math.Sin(2 * Math.PI / 3) / 3, 2)), 0.2, 1.5) |
                        
                        IneqTreeCylinder(0.75, 1.5 * Math.Sin(2 * Math.PI / 3), 0, 2.25 / 3, 1.5 * Math.Sin(2 * Math.PI / 3) / 3, Math.Sqrt(2.25 - Math.Pow(2 * 1.5 * Math.Sin(2 * Math.PI / 3) / 3, 2)), 0.2, 1.5)
            };

            samples["Helix"] = new IneqMesh
            {
                X0 = -2.0,
                Y0 = -2.0,
                Z0 = 0.0,
                X1 = 2.0,
                Y1 = 2.0,
                Z1 = 8.0,
                D = 0.15d,
                Boxed = true,
                IneqTree =
                        ((IneqTree)((x, y, z) => HelixDistance.Squared(x,y,z, 1, 0.2) - 0.35d))
            };

            //samples["* Bear"] = new Func<IneqMesh>(() => Bear());
            //samples["* Muse"] = new Func<IneqMesh>(() => Muse());

            //samples["* Bunny"] = new Func<IneqMesh>(() => Bunny());
            //samples["* BunnyHead"] = new Func<IneqMesh>(() => BunnyHead());
            //samples["* BunnyMuse"] = new Func<IneqMesh>(() => BunnyMuse());

            //samples["* Head"] = new Func<IneqMesh>(() => Head());
            //samples["* Head1"] = new Func<IneqMesh>(() => Head1());

            //samples["* Moai"] = new Func<IneqMesh>(() => Moai());
            //samples["* Moai1"] = new Func<IneqMesh>(() => Moai1());

            //samples["** Tooth"] = new Func<IneqMesh>(() => Tooth());

            //samples["** Skull and brain"] = new Func<IneqMesh>(() => Skull4());
            //samples["** Skull part"] = new Func<IneqMesh>(() => Skull(1));
        }

        public IneqMesh this[string name]
        {
            get
            {
                if (!(samples[name] is IneqMesh))
                {
                    //Cursor.Current = Cursors.WaitCursor;
                    samples[name] = ((Func<IneqMesh>)samples[name])();
                    //Cursor.Current = Cursors.Default;
                }
                return samples[name] as IneqMesh;
            }
        }


        private IneqTree IneqTreeBall(double x0, double y0, double z0, double r)
        {
            return new IneqTree((x, y, z) => (x - x0) * (x - x0) + (y - y0) * (y - y0) + (z - z0) * (z - z0) - r * r);
        }

        private IneqTree IneqTreeAnuloid(double x0, double y0, double z0, double r, double R, double nx, double ny, double nz, bool outside = false)
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

        private IneqTree IneqTreeCylinder(double x1, double y1, double z1, double x2, double y2, double z2, double r, double d)
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
                    double vp = (x - (x1 + x2) /2) * nx + (y - (y1 + y2) /2) * ny + (z - (z1 + z2) /2) * nz;
                    return Math.Abs(vp) - d;
                });
        }
        
        private IneqTree IneqTreeBalls(int count, double R, double r)
        {
            IneqTree res = new IneqTree((x, y, z) => 1);

            for (int i = 0; i < count; i++)
            {
                double x0 = R * Math.Cos(i * 2 * Math.PI / count);
                double y0 = R * Math.Sin(i * 2 * Math.PI / count);

                res = res | ((x, y, z) => (x - x0) * (x - x0) + (y - y0) * (y - y0) + z * z - r * r);
            }

            return res;
        }

        private IneqTree IneqTreeAnuloids(int count, double z0, double z1, double r0, double r1, double R0, double R1)
        {
            IneqTree res = new IneqTree((x, y, z) => 1);

            for (int i = 0; i < count; i++)
            {
                double z = z0 + i * (z1 - z0) / (count - 1);
                double r = r0 + i * (r1 - r0) / (count - 1);
                double R = R0 + i * (R1 - R0) / (count - 1);

                res = res | IneqTreeAnuloid(0, 0, z, r, R, 0, 0, 1);
            }

            return res;
        }

        private IneqTree IneqTreeAnuloids1(int count, double RR, double R,double r)
        {
            IneqTree res = new IneqTree((x, y, z) => 1);

            for (int i = 0; i < count; i++)
            {
                double rx = Math.Cos(i * 2 * Math.PI / count);
                double ry = Math.Sin(i * 2 * Math.PI / count);


                double x0 = RR* rx;
                double y0 = RR * ry;

                res = res | IneqTreeAnuloid(x0, y0, 0, r, R, -ry, rx, 0);
            }

            return res;
        }

        private IneqTree IneqTreeAntiAnuloids(int count, double RR, double R, double r)
        {
            IneqTree res = new IneqTree((x, y, z) => -1);

            for (int i = 0; i < count; i++)
            {
                double rx = Math.Cos(i * 2 * Math.PI / count);
                double ry = Math.Sin(i * 2 * Math.PI / count);


                double x0 = RR * rx;
                double y0 = RR * ry;

                res = res & IneqTreeAnuloid(x0, y0, 0, r, R, -ry, rx, 0, true);
            }

            return res;
        }

        private IneqTree IneqTreeCylinders(int count, double R, double r)
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

        private IneqTree IneqTreeCylinders1(int count, double R, double r)
        {
            IneqTree res = new IneqTree((x, y, z) => 1);

            for (int i = 0; i < count; i++)
            {
                double x0 = Math.Cos(i * 2 * Math.PI / count);
                double y0 = Math.Sin(i * 2 * Math.PI / count);

                res = res | (IneqTree) (((x, y, z) => 
                {
                    double vp = x * x0 + y * y0; 
                    return vp < 0 ? 1 : (x-vp*x0)*(x-vp*x0) + (y-vp*y0)*(y-vp*y0) + z*z - r*r;
                }) &
                (IneqTree)
                    ((x, y, z) => 
                    {
                        double vp = x * x0 + y * y0;
                        return -vp  + R;
                    }
                ));
            }

            return res;
        }

        private IneqTree IneqTreePlanes(int count)
        {
            IneqTree res = new IneqTree((x, y, z) => -1);

            for (int i = 0; i < count; i++)
            {
                double x0 = 0.75d * Math.Cos(i * 2 * Math.PI / count);
                double y0 = 0.75d * Math.Sin(i * 2 * Math.PI / count);
                double z0 = 0;

                res = res & ((x, y, z) => (x - x0) * x0 + (y - y0) * y0 + (z - z0) * z0 );
            }

            return res;
        }

        private IneqTree IneqTreePlanes1(int count)
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

        ////http://www.karlson.ru/csrbf/index.php?id=2
        //ComputeDistance computeDistanceBunny;
        //private IneqMesh Bunny()
        //{
        //    computeDistanceBunny = computeDistanceBunny ?? new ComputeDistance(@"VTKData\bunny.vtk", 30);

        //    return new IneqMesh
        //    {
        //        X0 = computeDistanceBunny.GetBounds()[0],
        //        Y0 = computeDistanceBunny.GetBounds()[2],
        //        Z0 = computeDistanceBunny.GetBounds()[4],
        //        X1 = computeDistanceBunny.GetBounds()[1],
        //        Y1 = computeDistanceBunny.GetBounds()[3],
        //        Z1 = computeDistanceBunny.GetBounds()[5],
        //        D = (computeDistanceBunny.GetBounds()[1] - computeDistanceBunny.GetBounds()[0]) / 50,
        //        Boxed = false,
        //        EvalTotalFunc = ((x, y, z) => x * x + (y - 0.9) * (y - 0.9) + z * z),
        //        IneqTree =
        //                ((IneqTree)((x, y, z) => computeDistanceBunny.GetVal(x, y, z)) & 
        //                ((x, y, z) => y - 0.96)) |
        //                (IneqTree)((x, y, z) => x * x + (y - 0.9) * (y - 0.9) + z * z - 0.02)
        //    };
        //}


        //ComputeDistance computeDistanceBunnyHead;
        //private IneqMesh BunnyHead()
        //{
        //    computeDistanceBunnyHead = computeDistanceBunnyHead ?? new ComputeDistance(@"VTKData\Bunny.vtk", 50, 0.1d, -1, false);

        //    return new IneqMesh
        //    {
        //        X0 = computeDistanceBunnyHead.GetBounds()[0],
        //        Y0 = computeDistanceBunnyHead.GetBounds()[2] + 0.3,
        //        Z0 = computeDistanceBunnyHead.GetBounds()[4],
        //        X1 = computeDistanceBunnyHead.GetBounds()[1] - 0.6,
        //        Y1 = computeDistanceBunnyHead.GetBounds()[3],
        //        Z1 = computeDistanceBunnyHead.GetBounds()[5],
        //        D = (computeDistanceBunnyHead.GetBounds()[1] - computeDistanceBunnyHead.GetBounds()[0]) / 80,
        //        Boxed = true,
        //        EvalTotalFunc = ((x, y, z) => x),
        //        IneqTree =
        //                ((IneqTree)((x, y, z) => computeDistanceBunnyHead.GetVal(x, y, z)))
        //    };
        //}

        //private IneqMesh BunnyMuse()
        //{
        //    computeDistanceBunny = computeDistanceBunny ?? new ComputeDistance(@"VTKData\bunny.vtk", 30); //!!!
        //    computeDistanceMuse = computeDistanceMuse ?? new ComputeDistance(@"VTKData\muse.3ds", 30);
        //    /*
        //        [0]: -0.62066001929342751
        //        [1]: 0.404340011253953
        //        [2]: 0.19936600588262082
        //        [3]: 1.2155999664217234
        //        [4]: -0.40989200510084628
        //        [5]: 0.3901499930769205
        //     */

        //    return new IneqMesh
        //    {
        //        X0 = computeDistanceBunny.GetBounds()[0],
        //        Y0 = computeDistanceBunny.GetBounds()[2],
        //        Z0 = computeDistanceBunny.GetBounds()[4],
        //        X1 = computeDistanceBunny.GetBounds()[1],
        //        Y1 = computeDistanceBunny.GetBounds()[3]+0.5,
        //        Z1 = computeDistanceBunny.GetBounds()[5],
        //        D = (computeDistanceBunny.GetBounds()[1] - computeDistanceBunny.GetBounds()[0]) / 40,
        //        Boxed = false,
        //        EvalTotalFunc = ((x, y, z) => x),
        //        IneqTree =
        //                ((IneqTree)((x, y, z) => computeDistanceBunny.GetVal(x, y, z)) &
        //                ((x, y, z) => y - 0.9)) |
        //                (IneqTree)((x, y, z) => computeDistanceMuse.GetVal(400 * (x + 0.45), -400 * (z - 0.1), 400 * (y - 1.08))) |
        //                //(IneqTree)((x, y, z) => (x + 0.4) * (x + 0.4) + (y - 1) * (y - 1) + (z-0.2) * (z-0.2) - 0.02)
        //                ((x, y, z) => (x + 0.48) * (x + 0.48) + 0.8*(y - 1.35) * (y - 1.35) + (z - 0.05) * (z - 0.05) - 0.003)                        
        //    };
        //}


        //ComputeDistance computeDistanceHead; 
        //private IneqMesh Head()
        //{
        //    computeDistanceHead = computeDistanceHead ?? new ComputeDistance(@"VTKData\head.vtk", 40, 0.2d);

        //    return new IneqMesh
        //    {
        //        X0 = computeDistanceHead.GetBounds()[0],
        //        Y0 = computeDistanceHead.GetBounds()[2],
        //        Z0 = computeDistanceHead.GetBounds()[4],
        //        X1 = computeDistanceHead.GetBounds()[1],
        //        Y1 = computeDistanceHead.GetBounds()[3],
        //        Z1 = computeDistanceHead.GetBounds()[5],
        //        D = (computeDistanceHead.GetBounds()[1] - computeDistanceHead.GetBounds()[0]) / 40,
        //        Boxed = false,
        //        EvalTotalFunc = ((x, y, z) => x),
        //        IneqTree =
        //                ((IneqTree)((x, y, z) => computeDistanceHead.GetVal(x, y, z)) &
        //                ((x, y, z) => -y + 2.5)) &
        //                ((x, y, z) => -x * x - z * z + 3)
        //    };
        //}

        //private IneqMesh Head1()
        //{
        //    computeDistanceHead = computeDistanceHead ?? new ComputeDistance(@"VTKData\head.vtk", 40, 0.2d);

        //    return new IneqMesh
        //    {
        //        X0 = computeDistanceHead.GetBounds()[0],
        //        Y0 = computeDistanceHead.GetBounds()[2],
        //        Z0 = computeDistanceHead.GetBounds()[4],
        //        X1 = computeDistanceHead.GetBounds()[1],
        //        Y1 = computeDistanceHead.GetBounds()[3],
        //        Z1 = computeDistanceHead.GetBounds()[5],
        //        D = (computeDistanceHead.GetBounds()[1] - computeDistanceHead.GetBounds()[0]) / 40,
        //        Boxed = false,
        //        EvalTotalFunc = ((x, y, z) => x),
        //        IneqTree =
        //                (((IneqTree)((x, y, z) => computeDistanceHead.GetVal(x, y, z)) &
        //                ((x, y, z) => -computeDistanceHead.GetVal(1.3 * x,  1.3 * (y - 3),  1.3 * z))
        //                ) &
        //                ((x, y, z) => -y + 3.5) &
        //                ((x, y, z) => x - 1)) 
        //    };
        //}


        //ComputeDistance computeDistanceMoai;
        //private IneqMesh Moai()
        //{
        //    computeDistanceMoai = computeDistanceMoai ?? new ComputeDistance(@"VTKData\moai.vtk", 20);
        //    return new IneqMesh
        //    {
        //        X0 = computeDistanceMoai.GetBounds()[0],
        //        Y0 = computeDistanceMoai.GetBounds()[2] - 0.2,
        //        Z0 = computeDistanceMoai.GetBounds()[4] - 0.3,
        //        X1 = computeDistanceMoai.GetBounds()[1],
        //        Y1 = computeDistanceMoai.GetBounds()[3],
        //        Z1 = computeDistanceMoai.GetBounds()[5] + 0.3,
        //        D = (computeDistanceMoai.GetBounds()[1] - computeDistanceMoai.GetBounds()[0]) / 20,
        //        Boxed = false,
        //        EvalTotalFunc = ((x, y, z) => x),
        //        IneqTree =
        //                ((IneqTree)((x, y, z) => computeDistanceMoai.GetVal(x, y, -(z + 0.25)))) |
        //                ((x, y, z) => computeDistanceMoai.GetVal(x, y, z - 0.25)) |
        //                (
        //                    (IneqTree)((x, y, z) => Math.Abs(x) - 0.25) &
        //                    ((x, y, z) => Math.Abs(z) - 0.5) &
        //                    ((x, y, z) => Math.Abs(y + 0.51) - 0.1)                        
        //                )
        //    };
        //}

        //private IneqMesh Moai1()
        //{
        //    computeDistanceMoai = computeDistanceMoai ?? new ComputeDistance(@"VTKData\moai.vtk", 20);
        //    return new IneqMesh
        //    {
        //        X0 = computeDistanceMoai.GetBounds()[0] - 0.3,
        //        Y0 = computeDistanceMoai.GetBounds()[2] - 0.2,
        //        Z0 = computeDistanceMoai.GetBounds()[4] - 0.3,
        //        X1 = computeDistanceMoai.GetBounds()[1] + 0.3,
        //        Y1 = computeDistanceMoai.GetBounds()[3],
        //        Z1 = computeDistanceMoai.GetBounds()[5] + 0.3,
        //        D = (computeDistanceMoai.GetBounds()[1] - computeDistanceMoai.GetBounds()[0]) / 20,
        //        Boxed = false,
        //        EvalTotalFunc = ((x, y, z) => x),
        //        IneqTree =
        //                ((IneqTree)((x, y, z) => computeDistanceMoai.GetVal(x, y, -(z + 0.25)))) |
        //                ((x, y, z) => computeDistanceMoai.GetVal(x, y, z - 0.25)) |
        //                ((x, y, z) => computeDistanceMoai.GetVal(-z , y, (x - 0.25))) |
        //                ((x, y, z) => computeDistanceMoai.GetVal(-z, y, (-(x + 0.25)))) |
        //                (
        //                    (IneqTree)((x, y, z) => Math.Abs(x) - 0.5) &
        //                    ((x, y, z) => Math.Abs(z) - 0.5) &
        //                    ((x, y, z) => Math.Abs(y + 0.51) - 0.1)
        //                )
        //    };
        //}

        ////http://archive3d.net
        //ComputeDistance[] computeDistanceSkullParts = new ComputeDistance[5];
        //private IneqMesh Skull(int part)
        //{
        //    computeDistanceSkullParts[part] = computeDistanceSkullParts[part] ?? new ComputeDistance(@"VTKData\skull.3DS", 60, 0.25d, part, true);

        //    ComputeDistance computeDistanceSkull = computeDistanceSkullParts[part];
        //    return new IneqMesh
        //    {
        //        X0 = computeDistanceSkull.GetBounds()[0],
        //        Y0 = computeDistanceSkull.GetBounds()[2],
        //        Z0 = computeDistanceSkull.GetBounds()[4],
        //        X1 = computeDistanceSkull.GetBounds()[1],
        //        Y1 = computeDistanceSkull.GetBounds()[3],
        //        Z1 = computeDistanceSkull.GetBounds()[5],
        //        D = (computeDistanceSkull.GetBounds()[1] - computeDistanceSkull.GetBounds()[0]) / 60,
        //        Boxed = false,
        //        EvalTotalFunc = ((x, y, z) => x),
        //        IneqTree =
        //                ((IneqTree)((x, y, z) => computeDistanceSkull.GetVal(x, y, z)))
        //    };
        //}

        //private IneqMesh Skull4()
        //{
        //    computeDistanceSkullParts[4] = computeDistanceSkullParts[4] ?? new ComputeDistance(@"VTKData\skull.3DS", 60, 0.25d, 4, true);

        //    ComputeDistance computeDistanceSkull = computeDistanceSkullParts[4];
        //    return new IneqMesh
        //    {
        //        X0 = computeDistanceSkull.GetBounds()[0],
        //        Y0 = computeDistanceSkull.GetBounds()[2],
        //        Z0 = computeDistanceSkull.GetBounds()[4],
        //        X1 = computeDistanceSkull.GetBounds()[1],
        //        Y1 = computeDistanceSkull.GetBounds()[3],
        //        Z1 = computeDistanceSkull.GetBounds()[5],
        //        D = (computeDistanceSkull.GetBounds()[1] - computeDistanceSkull.GetBounds()[0]) / 60,
        //        Boxed = false,
        //        EvalTotalFunc = ((x, y, z) => x),
        //        IneqTree =
        //                ((IneqTree)((x, y, z) => computeDistanceSkull.GetVal(x, y, z)) &
        //                    ((x, y, z) => -x * x - y * y + 25) &
        //                    ((x, y, z) => -(x + 1) * (x + 1) - y * y - (z - 3) * (z - 3)  + 50)
        //                ) |
        //            (
        //                (IneqTree)
        //                    ((x, y, z) => (x + 1) * (x + 1) + y * y + (z - 3) * (z - 3) - (5 * (Math.Cos(2 * x) + Math.Cos(2 * y) + Math.Cos(2 * z)) + 38)) |
        //                    ((x, y, z) => (x + 1) * (x + 1) + y * y + (z - 3) * (z - 3) - (5 * (Math.Sin(2 * x) + Math.Sin(2 * y) + Math.Sin(2 * z)) + 38))
        //            )
        //    };
        //}


        //ComputeDistance computeDistanceMuse;
        //private IneqMesh Muse()
        //{
        //    /* Muse:
        //        [0]: -50.634556484222415
        //        [1]: 50.634556484222415
        //        [2]: -45.083035182952884
        //        [3]: 45.083035182952884
        //        [4]: -113.49074621200562
        //        [5]: 113.49075384140015             
        //     */

        //    computeDistanceMuse = computeDistanceMuse ?? new ComputeDistance(@"VTKData\muse.3ds", 30);
        //    return new IneqMesh
        //    {
        //        X0 = computeDistanceMuse.GetBounds()[0],
        //        Y0 = computeDistanceMuse.GetBounds()[2] - 5,
        //        Z0 = computeDistanceMuse.GetBounds()[4] - 20,
        //        X1 = computeDistanceMuse.GetBounds()[1],
        //        Y1 = computeDistanceMuse.GetBounds()[3] + 5,
        //        Z1 = computeDistanceMuse.GetBounds()[5] + 30,
        //        D = (computeDistanceMuse.GetBounds()[1] - computeDistanceMuse.GetBounds()[0]) / 30,
        //        Boxed = false,
        //        EvalTotalFunc = ((x, y, z) => x * x + 0.30 * (z + 90) * (z + 90) + y * y),
        //        IneqTree =
        //                ((IneqTree)((x, y, z) => computeDistanceMuse.GetVal(x, y, z)) |
        //                    ((x, y, z) => (x+10)*(x+10) +1.2 * (y-10)*(y-10) + 0.8 * (z-110)*(z-110) -600) |
        //                    ((IneqTree)((x, y, z) => 0.8 * x * x + 0.8 * (z + 90) * (z + 90) - 300) &
        //                    ((x, y, z) => Math.Abs(y) - 45))
        //                )
        //    };
        //}


        //ComputeDistance computeDistanceBear;
        //private IneqMesh Bear()
        //{
        //    computeDistanceBear = computeDistanceBear ?? new ComputeDistance(@"VTKData\Bear.3ds", 30); //30
        //    return new IneqMesh
        //    {
        //        X0 = computeDistanceBear.GetBounds()[0],
        //        Y0 = computeDistanceBear.GetBounds()[2],
        //        Z0 = computeDistanceBear.GetBounds()[4],
        //        X1 = computeDistanceBear.GetBounds()[1],
        //        Y1 = computeDistanceBear.GetBounds()[3] + 10,
        //        Z1 = computeDistanceBear.GetBounds()[5],
        //        D = (computeDistanceBear.GetBounds()[1] - computeDistanceBear.GetBounds()[0]) / 40, //40
        //        Boxed = false,
        //        EvalTotalFunc = ((x, y, z) => /*x*/ (x + 6) * (x + 6) + (y - 2.7) * (y - 2.7) + (z - 7.6) * (z - 7.6)),
        //        IneqTree =
        //                (
        //                    (IneqTree)((x, y, z) => computeDistanceBear.GetVal(x, y, z)) |
        //                    ((x, y, z) => (x + 6) * (x + 6) + (y - 2.7) * (y - 2.7) + (z - 7.6) * (z - 7.6) - 16) |
        //                    (((x, y, z) =>
        //                        {
        //                            y -= 25;
        //                            x -= 15;
        //                            z += 15;

        //                            double n = Math.Sqrt(2.0d) / 2.0d;
        //                            double p =  n * x + n * z;

        //                            return (x - p * n) * (x - p * n) + (z - p * n) * (z - p * n) + y * y - 25;
        //                        }) &
        //                        (IneqTree)((x, y, z) => Math.Abs(x + z)-20))
        //                )
        //    };
        //}


        //ComputeDistance computeDistanceTooth;
        //private IneqMesh Tooth()
        //{
        //    computeDistanceTooth = computeDistanceTooth ?? new ComputeDistance(@"VTKData\tooth.3ds", 30);
        //    /* Tooth:
        //        [0]: -42.248913860321046
        //        [1]: 38.6016131401062
        //        [2]: -40.265145397186281
        //        [3]: 39.224057292938234
        //        [4]: -1.5288268089294434
        //        [5]: 57.3336688041687             
        //     */
        //    return new IneqMesh
        //    {
        //        X0 = computeDistanceTooth.GetBounds()[0],
        //        Y0 = computeDistanceTooth.GetBounds()[2],
        //        Z0 = computeDistanceTooth.GetBounds()[4],
        //        X1 = computeDistanceTooth.GetBounds()[1],
        //        Y1 = computeDistanceTooth.GetBounds()[3],
        //        Z1 = computeDistanceTooth.GetBounds()[5],
        //        D = (computeDistanceTooth.GetBounds()[1] - computeDistanceTooth.GetBounds()[0]) / 40,
        //        Boxed = false,
        //        EvalTotalFunc = ((x, y, z) => x * x + y * y ),
        //        IneqTree =
        //                ((IneqTree)((x, y, z) => computeDistanceTooth.GetVal(x, y, z))) &
        //                ((x, y, z) => -(x + 15) * (x + 15) - (y - 15) * (y - 15) + 25) &
        //                ((x, y, z) => -(x - 4) * (x - 4) - (y - 29) * (y - 29) + 49) &
        //                ((x, y, z) => -(x + 5) * (x + 5) - (y + 25) * (y + 25) + 25)
        //    };
        //}

        //if (Mesh == "Cone 1")
        //{
        //    Point p = new Point(0.12, -0.171, 0.037);
        //    Point np = null;
        //    double d = Double.MaxValue;
        //    foreach (var mp in ineqMesh.Points.Where(pp => pp.BoundaryCount == 1))
        //    {
        //        if (mp.Distance(p) < d)
        //        {
        //            d = mp.Distance(p);
        //            np = mp;
        //        }
        //    }

        //    np.MoveTo(p, false);
        //    np.Boundary[np.Boundary.Length - 1] = true;
        //    ineqMesh.Jiggle(3);

        //    ineqMesh.RefineBoundaryTriangles(np.Tetrahedrons.SelectMany(t => t.Triangles()).Where(t => t.All(pp => pp.BoundaryCount > 0)));
        //    ineqMesh.Jiggle(3);
        //    ineqMesh.RefineBoundaryTriangles(ineqMesh.Triangles.Where(tt => tt.All(pp => pp.Boundary[0])));
        //    ineqMesh.Jiggle(3);
        //}

        //ineqMesh.RefineBoundaryTriangles(ineqMesh.Triangles.Where(tt => tt.All(p => p.Boundary[5])));
        //ineqMesh.Jiggle(3);
    }
}