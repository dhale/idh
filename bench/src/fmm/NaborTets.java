package fmm;

/**
 * Computes indices of tets for each neighbor sample in a 26-neighborhood.
 * This class only prints a table that is needed elsewhere.
 */
public class NaborTets {

  public static void main(String[] args) {
    System.out.println("{");
    for (int k3=-1; k3<=1; ++k3) {
      for (int k2=-1; k2<=1; ++k2) {
        for (int k1=-1; k1<=1; ++k1) {
          if (k1==0 && k2==0 && k3==0) continue;
          System.out.print("  {");
          for (int kt=0; kt<48; ++kt) {
            int k11 = K11[kt];
            int k12 = K12[kt];
            int k13 = K13[kt];
            int k21 = K21[kt];
            int k22 = K22[kt];
            int k23 = K23[kt];
            int k31 = K31[kt];
            int k32 = K32[kt];
            int k33 = K33[kt];
            if (k11==k1 && k12==k2 && k13==k3 ||
                k21==k1 && k22==k2 && k23==k3 ||
                k31==k1 && k32==k2 && k33==k3)
              System.out.printf("%2d,",kt);
          }
          System.out.println("},");
        }
      }
    }
    System.out.println("};");
  }

  // Sample index offsets for vertices X1 of the 48 nabor tetrahedra.
  private static final int[] K11 = { 1, 1, 0,-1,-1,-1, 0, 1,
                                     1, 1, 0,-1,-1,-1, 0, 1,
                                     1, 1, 0,-1,-1,-1, 0, 1,
                                     1, 1, 0,-1,-1,-1, 0, 1,
                                    -1,-1,-1,-1,-1,-1,-1,-1,
                                     1, 1, 1, 1, 1, 1, 1, 1};
  private static final int[] K12 = { 0, 1, 1, 1, 0,-1,-1,-1,
                                     0, 1, 1, 1, 0,-1,-1,-1,
                                    -1,-1,-1,-1,-1,-1,-1,-1,
                                     1, 1, 1, 1, 1, 1, 1, 1,
                                     1, 1, 0,-1,-1,-1, 0, 1,
                                     1, 1, 0,-1,-1,-1, 0, 1};
  private static final int[] K13 = {-1,-1,-1,-1,-1,-1,-1,-1,
                                     1, 1, 1, 1, 1, 1, 1, 1,
                                     0, 1, 1, 1, 0,-1,-1,-1,
                                     0, 1, 1, 1, 0,-1,-1,-1,
                                     0, 1, 1, 1, 0,-1,-1,-1,
                                     0, 1, 1, 1, 0,-1,-1,-1};

  // Sample index offsets for vertices X2 of the 48 nabor tetrahedra.
  private static final int[] K21 = { 1, 0,-1,-1,-1, 0, 1, 1,
                                     1, 0,-1,-1,-1, 0, 1, 1,
                                     1, 0,-1,-1,-1, 0, 1, 1,
                                     1, 0,-1,-1,-1, 0, 1, 1,
                                    -1,-1,-1,-1,-1,-1,-1,-1,
                                     1, 1, 1, 1, 1, 1, 1, 1};
  private static final int[] K22 = { 1, 1, 1, 0,-1,-1,-1, 0,
                                     1, 1, 1, 0,-1,-1,-1, 0,
                                    -1,-1,-1,-1,-1,-1,-1,-1,
                                     1, 1, 1, 1, 1, 1, 1, 1,
                                     1, 0,-1,-1,-1, 0, 1, 1,
                                     1, 0,-1,-1,-1, 0, 1, 1};
  private static final int[] K23 = {-1,-1,-1,-1,-1,-1,-1,-1,
                                     1, 1, 1, 1, 1, 1, 1, 1,
                                     1, 1, 1, 0,-1,-1,-1, 0,
                                     1, 1, 1, 0,-1,-1,-1, 0,
                                     1, 1, 1, 0,-1,-1,-1, 0,
                                     1, 1, 1, 0,-1,-1,-1, 0};

  // Sample index offsets for vertices X3 of the 48 nabor tetrahedra.
  private static final int[] K31 = { 0, 0, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0, 0,
                                    -1,-1,-1,-1,-1,-1,-1,-1,
                                     1, 1, 1, 1, 1, 1, 1, 1};
  private static final int[] K32 = { 0, 0, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0, 0,
                                    -1,-1,-1,-1,-1,-1,-1,-1,
                                     1, 1, 1, 1, 1, 1, 1, 1,
                                     0, 0, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0, 0};
  private static final int[] K33 = {-1,-1,-1,-1,-1,-1,-1,-1,
                                     1, 1, 1, 1, 1, 1, 1, 1,
                                     0, 0, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0, 0};
}

