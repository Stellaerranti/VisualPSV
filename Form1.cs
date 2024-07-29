using System;
using System.Collections;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace VisualPSV
{
    public partial class Form1 : Form
    {
        //Main lists
        private List<double> Dec = new List<double>();
        private List<double> Inc = new List<double>();
        private List<double> R = new List<double>();

        private List<double> X = new List<double>();
        private List<double> Y = new List<double>();
        private List<double> Z = new List<double>();

        public Form1()
        {
            InitializeComponent();
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            Dec.Clear();
            Inc.Clear();
            R.Clear();

            X.Clear();
            Y.Clear();
            Z.Clear();

            TK03Session();

        }

        //Convert directionsl coords to cartesizn
        private List<double> dir2cart(List<double> dir) 
        {
            double rad = Math.PI / 180;

            double dec = dir[0];
            double inc = dir[1];
            double r = dir[2];

            double x = r * Math.Cos(dec) * Math.Cos(inc);
            double y = r * Math.Sin(dec) * Math.Cos(inc);
            double z = r * Math.Sin(inc);

            return dir;
        }

        //Convert cartesizn coords to directional
        private List<double> cart2dir(List<double> cart)
        {
            double x = cart[0];
            double y = cart[1];
            double z = cart[2];

            double rad = Math.PI / 180;

            double r = Math.Sqrt(Math.Pow(x,2)+ Math.Pow(y, 2)+ Math.Pow(z, 2));

            double dec = (Math.Atan2(y, x) / rad) % 360;

            double inc = Math.Asin(z/r)/rad;

            List<double> dir = new List<double>() { dec, inc, r };

            return dir;
        }

        //Get normmal distibution using a  Box - Muller transform
        private double RandomNormal(double mean, double std) 
        {
            Random rand = new Random(); //reuse this if you are generating many
            double u1 = 1.0 - rand.NextDouble(); //uniform(0,1] random doubles
            double u2 = 1.0 - rand.NextDouble();
            double randStdNormal = Math.Sqrt(-2.0 * Math.Log(u1)) *
                         Math.Sin(2.0 * Math.PI * u2); //random normal(0,1)
            double randNormal =
                         mean + std * randStdNormal; //random normal(mean,std^2)

            return randNormal;
        }

        private void TK03Session()
        {
            Dec.Clear();
            Inc.Clear();
            R.Clear();

            X.Clear();
            Y.Clear();
            Z.Clear();

            List<double> gh = new List<double>();
            List<double> sv = new List<double>();

            List<double> vecCoords = new List<double>();

            double lat = 60;
            double lon = 60;

            double colat = 0;



            int k = 0;
            int n = 20;

            while (k<n)
            {
                gh = getCoeffsTK03(8, 18e3,0,0, 3.8, 2.4);
                
                while (gh.Count < 120)
                {gh.Add(0);}

                for (int i =0; i<gh.Count;i++)
                { sv.Add(0); }

                colat = 90 - lat;

                vecCoords = getCoords(gh, sv, 0, 0, 1, 0, colat, lon);

                File.AppendAllText("test.txt", vecCoords[0].ToString("F4") +" " + vecCoords[1].ToString("F4") + " " + vecCoords[2].ToString("F4") + " " + vecCoords[3].ToString("F4")+ "\n");
                k++;
            }



            /*
             
                    using (StreamWriter outputFile = new StreamWriter(Path.Combine(docPath, "test.txt")))
        {
            foreach (string line in lines)
                outputFile.WriteLine(line);
        }
             for (int row = 0; row < 13; row++)
{
    for (int col = 0; col < 40; col++)
    {
        outfile.Write("{0},", monthData[row, col]); // note you had missed a comma
        outfile.Write("{0}", Environment.Newline);
    }
}
             */
        }

        //Returns sigma, corresponding to hamonic l
        private double SL(int l, double alpha = 2.7)
        {

            double a2 = Math.Pow(alpha, 2);

            double ca = 0.547;

            double s_l = Math.Sqrt(a2*Math.Pow(ca,2*l)/((l+1)*(2*l+1)));
            return s_l;
        }

        // This functions returns h and g coeffs for TK03
        private List<double> getCoeffsTK03(int terms, double G1, double G2, double G3, double beta, double afact)
        {
            List<double> gh = new List<double>();

            double alpha = G1/afact;

            double s1 = SL(1, alpha);
            
            double s10 = s1 * beta;

            gh.Add(RandomNormal(G1,s10));
            gh.Add(RandomNormal(0,s1));
            gh.Add(RandomNormal(0, s1));

            double s = 0;
            double o = 0;

            for (int l = 2; l < (terms+1);l++)
            {
                for (int m = 0; m < l+1; m++)
                {
                    o = 0;

                    if(l == 2 && m == 0)
                    {o = G2;}

                    if (l == 3 && m == 0)
                    { o = G3; }

                    s = SL(l,alpha);

                    if ((l - m) % 2 == 1)
                    { s = s * beta; }

                    gh.Add(RandomNormal(o,s));

                    if(m != 0)
                    {gh.Add(RandomNormal(0, s));}
                }
            }

            return gh;
        }

        //Get coords from gh-s and sv-s
        private List<double> getCoords(List<double> gh, List<double> sv, double b, double date, int itype, double altitude, double colatitude, double longitude)
        {
            double x,y,z,f;
            x=y=z=f= 0;

            double[] p = new double[66];
            double[] q = new double[66];
            double[] sl = new double[10];
            double[] cl = new double[10];

            //For redundancy
            Array.Clear(p, 0, p.Length);
            Array.Clear(q, 0, q.Length);
            Array.Clear(sl, 0, sl.Length);
            Array.Clear(cl, 0, cl.Length);

            double t = date - b;

            double one, two, three;
            one = two = three = 0;

            one = colatitude * 0.0174532925;

            double ct = Math.Cos(one);
            double st = Math.Sin(one);

            cl[0] = ct;
            sl[0] = st;

            double sd = 0;
            double cd = 1;

            double r = altitude;
            double ratio = 0;
            double rr = 0;

            int l = 1;
            int ll = 0;
            int m = 1;
            int n = 0;
            int j = 0;
            int i = 0;

            double fn = 0;
            double fm = 0;
            double gm = 0;
            double gn = 0;
            
            

            double a2 = 40680925.0;
            double b2 = 40408585.0;

            double rho = 0;

            //transform geodectic to geocentric 
            if (itype == 1)
            {
                one = a2 * st * st;
                two = b2 * ct * ct;

                three = one + two;

                
                rho = Math.Sqrt(three);

                r = Math.Sqrt(r * (r+2*rho) + ((a2*one+b2*two)/three));

                cd = (altitude + rho) / r;
                sd = (a2 - b2)/ ct * rho *st/r;

                one = ct;
                ct = ct * cd - st * sd;
                st = st * cd + one * sd;
            }
            ratio = 6371.2 / r;
            rr = ratio*ratio;

            //Schmidt quasi-normal coefficients p and x(=q)

            p[0] = 1;
            p[2] = st;
            q[0] = 0.0;
            q[2] = ct;

            //Main cycle

            for(int k = 1; k<66; k++)
            {
                if(n < m)
                {
                    m = 0;
                    n++;
                    rr = rr * ratio;
                    fn = n;
                    gn = n - 1;
                }
                fm = m;

                if(k != 2)
                {
                    if (m == n)
                    {
                        one = Math.Sqrt(1 - 0.5/fm);
                        j = k - n - 1;
                        p[k] = one * st * p[j];
                        q[k] = one * (st * q[j]+ ct * p[j]);

                        cl[m - 1] = cl[m - 2] * cl[0] - sl[m - 2] * sl[0];
                        cl[m - 1] = sl[m - 2] * cl[0] + cl[m - 2] * sl[0];
                    }
                    else 
                    {
                        gm = m * m;

                        one = Math.Sqrt(fn * fn - gm);
                        two = Math.Sqrt(gn * gn - gm)/one;
                        three = (fn + gn) / one;

                        i = k - n;
                        j = i - n + 1;

                        p[k] = three * ct * p[i] - two * p[j];
                        q[k] = three * (ct * q[i] - st * p[i]) - two * q[j];
                    }
                }

                one = (gh[l - 1] + sv[ll + l - 1] * t) * rr;

                if(m != 0)
                {
                    //two = (gh[l] + sv[ll + l] * t) * rr;
                    two = gh[l] * rr;
                    three = one * cl[m - 1] + two * sl[m - 1];

                    x = x + three * q[k];
                    z = z - (fn + 1) * three * p[k];

                    if(st != 0)
                    {
                        y = y + (one * sl[m - 1] - two * cl[m - 1]) * fm * p[k]/st;
                    }
                    else
                    {
                        y = y + (one * sl[m - 1] - two * cl[m - 1]) * q[k] * ct;
                    }

                    l = l + 2;
                }
                else 
                {
                    x = x + one * q[k];
                    z = z - (fn + 1) * one * p[k];
                    l = l + 1;
                }

                m = m + 1;
            }

            one = x;

            x = x * cd + z * sd;
            z = z * cd - one * sd;
            f = Math.Sqrt(x * x + y * y + z * z);

            List<double> coords = new List<double>() {x,y,z,f};
            return coords;
        }


    }
}
