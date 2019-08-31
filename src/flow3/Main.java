package ksm;

import static java.lang.System.out;
import java.util.HashMap;

class Main
{
    private static String key(String s1, String s2) {
        return s1 + ":" + s2;
    }

    public static void main(String[] argv)
    {
        String[] d = {"d1", "d2", "d3", "d4"};
        int m = d.length;

        String[] o = {"o1", "o2", "o3"};
        int n = o.length;

        // pair-wise distances matrix
        HashMap<String, Integer> T = new HashMap();
        for (int i=0; i < n; i++)
            for (int j=0; j < n; j++)
                T.put(key(o[i], o[j]), i+j);

        for (int i=0; i < m; i++)
            for (int j=0; j < m; j++)
                T.put(key(d[i], d[j]), i+j);

        for (int i=0; i < n; i++)
            for (int j=0; j < m; j++)
                T.put(key(o[i], d[j]), i+j);

        for (int i=0; i < m; i++)
            for (int j=0; j < n; j++)
                T.put(key(d[i], o[j]), i+j);


        // initialize destination-subsequence length matrix l[][]
        // of size (n+1)x(n+1)
        int[][] l = new int[n+1][n+1];
        for (int xk=0; xk <= n; xk++) {
            for (int xk_1=0; xk_1 <= xk; xk_1++) {
                l[xk_1][xk] = 0;
                for (int i=xk_1+1; i <= xk-1; i++)
                    l[xk_1][xk] += T.get(key(d[i], d[i+1]));
            }
        }


    }
}
