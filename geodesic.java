import java.util.*;
import java.util.Arrays.*;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.ode.*;
import org.apache.commons.math3.ode.nonstiff.*;
import org.apache.commons.math3.ode.sampling.*;
import java.io.FileWriter;
import java.io.IOException;
import java.io.File;

class geodesic{
    //main class
    /*
    Main constants
    a and b are the Boyer-Lindquist coordinates and the Schwarzschild radius
    respectively
    */
    public static final double T = 1.0; public static final double R = 22.995;
    public static final double a = 1.0; public static final double b = 5.0;
    public static double[][] results;
    public static double[][] append(double[][] mat,double[] row){
        //Appends array to matrix
        double[][] newMatrix = new double[mat.length+1][mat[0].length];
        for(int j=0; j<mat[0].length; j++){
            for(int i=0; i<mat.length; i++){
                newMatrix[i][j] = mat[i][j];
            }
            newMatrix[mat.length][j]=row[j];
        }
        return newMatrix;
    }
    public static double[] linspace(double A, double B,int N){
        //Creates N evenly spaced numbers in [A,B]
        double M = (B-A)/N;
        double[] r = new double[N];
        for(int i = 0; i<N; i++){
            r[i] = A + i*M;
        }
        return r;
    }
    public static double[][] tphiMatrix(double r){
        /*Coefficient matrix for the linear relation
        of derivatives of time-angle*/
        double[][] R = {{-(1-b/r),-b*a/r},{-b*a/r, r*r+a*a+b*a*a/r}};
        /*for (int i = 0; i < R.length; i++) {
          for (int j = 0; j < R[i].length; j++) {
            System.out.print(R[i][j] + " ");
          }
          System.out.println();
        }*/
        System.out.println(r);
        return R;
    }
    public static double Delta(double r){
        //Delta function in the Kerr metric
        double R = a*a + r*r - r*b;
        return R;
    }
    public static RealVector tphi(double r){
        /*
        Taking solutions for the derivative of time and angle at a given value
        of the radius
        */
        RealMatrix coefficients = new Array2DRowRealMatrix(tphiMatrix(r),false);
        DecompositionSolver solver = new LUDecomposition(coefficients).getSolver();
        RealVector constants = new ArrayRealVector(new double[] {T,R}, false);
        RealVector solution = solver.solve(constants);
        return solution;
    }
    private static class Kerr implements FirstOrderDifferentialEquations {
        /*
        Equations of motion
        */
        public double RDP(double r){
          return (1-b/r)*R - a*b*T/r;
        }

        public double TDT(double r){
          return -(r*r+a*a+a*a*b/r)*T-a*b*R/r;
        }

        public double rDotSq(double r, double U, double V){
          return (1/r/r)*(-Delta(r)-T*U-R*V);
        }

        public int getDimension() {
            return 4;
        }
        public void computeDerivatives(double tau, double[] y, double[] yDot) {
            System.out.println(Arrays.toString(y));
            double r = y[2]; double t = y[0]; double phi = y[1];
            double D = Delta(r);
            double tt = TDT(r); double pp = RDP(r);
            //RealVector vec = tphi(r);
            //yDot[0] = vec.getEntry(0);
            //yDot[1] = vec.getEntry(1);
            yDot[0] = tt/D;
            yDot[1] = pp/D;
            yDot[2] = y[3];
            yDot[3] = (-1.0/2.0/r)*(rDotSq(r,tt/D,pp/D)+(1.0/r)*(2.0*r*(1.0-T*T)-b+b*(a*T+R)*(a*T+R)/r/r));
            //yDot[3] = (1.0/(2.0*r*r))*(-b*(yDot[0]-a*yDot[1])*(yDot[0]-a*yDot[1])+2*r*yDot[1]*yDot[1]-2*(D/r - (2*r+b))*(T*R-b*a/r + yDot[0]*yDot[1]*D));
            //yDot[3] = -(1/(D*2.0*r*r))*(b*r*r*yDot[2]*yDot[2]-2*r*r*r*yDot[2]*yDot[2] + D*D*b*(yDot[0]-a*yDot[1])*(yDot[0]-a*yDot[1])/r/r + r*(D*2*yDot[2]*yDot[2] - D*D*2*yDot[1]*yDot[1]));
            //yDot[3] = -rDotSq(r,tt,pp)/(r*r*r) + (1/(2*r*r))*(-2*a*b*T*R/r/r + b - 2*r + b*R*R/r/r + T*T*(2*r-a*a*b/r/r));
        }

    }

    public static void saveToCsv(){
        /*
        Used to export data to Python as csv for the plots
        */
        String[] order = new String[] {"t","phi","r","u","tau","x","y"};
        try{
           File file = new File("data.csv"); file.createNewFile();
           FileWriter writer = new FileWriter(file);
           for(String element : order){writer.write(","+element);}
           writer.write("\n");
           for(int j=0;j<results.length;j++){
               for (int i=0;i<results[0].length; i++){
                   writer.write(String.format(",%6.3e",results[j][i]));
               }
               writer.write("\n");
           }
           writer.close();
        }
        catch(IOException e){
           System.out.println(e.getMessage());
        }
   }

    public static void run(double tau0, double tauf, double steps, double r0){
        /*
        Solving differetial equation with DormandPrince853Integrator
        */
        FirstOrderIntegrator dp853 = new DormandPrince853Integrator(1.0e-10, 100.0, 1.0e-10, 1.0e-10);
        FirstOrderDifferentialEquations ode = new Kerr();
        double[] initialConditions = new double[] {0.0,0.0,r0,0.0,tau0,0.0,0.0};
        double realR0 = Math.sqrt(initialConditions[2]*initialConditions[2]+a*a);
        initialConditions[initialConditions.length-1]=realR0*Math.cos(initialConditions[1]);
        initialConditions[initialConditions.length-2]=realR0*Math.sin(initialConditions[1]);
        double[] y = new double[initialConditions.length-3]; for(int k = 0; k<initialConditions.length-3; k++){y[k]=initialConditions[k];} // initial state t,phi,r,rdot
        results = new double[][] {initialConditions};

        double deltaTau = (tauf-tau0)/steps;
        double[] out = new double[y.length+3];
        try{for(int i = 0; i<steps-1; i++){
            double tauI = tau0+i*deltaTau;
            double tauF = tauI + deltaTau;
            dp853.integrate(ode, tauI, y, tauF, y);
            for(int j=0;j<y.length;j++){
                out[j]=y[j];
            }
            out[y.length]=tauI;
            double realR = Math.sqrt(out[2]*out[2]+a*a);
            out[y.length+1]=realR*Math.cos(out[1]);
            out[y.length+2]=realR*Math.sin(out[1]);
            results = append(results, out);
        }}
        catch(Exception ignore){}
    }
    public static void main(String[] args){
       run(0.0,30000.0,1000.0, 100.0); //Initial conditions
       saveToCsv();
    }
}
