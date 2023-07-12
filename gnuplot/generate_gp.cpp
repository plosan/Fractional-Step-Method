#include <iostream>
#include <fstream>
#include <string>


using namespace std;

int main() {

    for(int i = 206; i <= 206; i += 5) {

        string title = "set multiplot layout 2,2 rowsfirst scale 1,1 title 'Re = 1000 (t = " + to_string(i) + " s, uref = 0.1 m/s, 129 x 129)'";
        string plotu = "plot 'data/u_129_129_1000_" + to_string(i) + ".txt' with image";
        string plotv = "plot 'data/v_129_129_1000_" + to_string(i) + ".txt' with image";
        string plotvel = "plot 'data/vel_129_129_1000_" + to_string(i) + ".txt' with image, 'data/vel_129_129_1000_" + to_string(i) + ".txt' using 1:2:($4/(20*sqrt(($4)**2 + ($5)**2))):($5/(20*sqrt(($4)**2 + ($5)**2))) every 5:5 with vectors lc -1 filled";
        string plotp = "plot 'data/p_129_129_1000_" + to_string(i) + ".txt' with image";

        string com[] = {
            "set terminal wxt size 1200,900",
            "\n",
            "\n",
            "set palette rgbformulae 22,13,-31",
            "unset key",
            "\n",
            title,
            "# --- GRAPH a",
            "\n",
            "set title 'U'",
            "\n",
            "set xrange [0:1]",
            "set yrange [0:1]",
            "set size ratio 1",
            "\n",
            "set xtics axis",
            "set xtics format '%.1f'",
            "set xlabel ('x [m]')",
            "\n",
            "\n",
            "set ytics axis",
            "set ytics format '%.1f'",
            "set ylabel ('y [m]')",
            "\n",
            "set cblabel ('u [m/s]')",
            "\n",
            "\n",
            plotu,
            "\n",
            "\n",
            "\n",
            "# --- GRAPH b",
            "\n",
            "set title 'v'",
            "\n",
            "set xrange [0:1]",
            "set yrange [0:1]",
            "set size ratio 1",
            "\n",
            "set xtics axis",
            "set xtics format '%.1f'",
            "set xlabel ('x [m]')",
            "\n",
            "set ytics axis",
            "set ytics format '%.1f'",
            "set ylabel ('y [m]')",
            "\n",
            "set cblabel ('v [m/s]')",
            "\n",
            "\n",
            "\n",
            plotv,
            "\n",
            "# --- GRAPH c",
            "\n",
            "\n",
            "set title 'velocity field'",
            "\n",
            "set xrange [0:1]",
            "set yrange [0:1]",
            "set size ratio 1",
            "\n",
            "set xtics axis",
            "set xtics format '%.1f'",
            "set xlabel ('x [m]')",
            "\n",
            "set ytics axis",
            "set ytics format '%.1f'",
            "set ylabel ('y [m]')",
            "\n",
            "set cblabel ('vel [m/s]')",
            "\n",
            "\n",
            "\n",
            plotvel,
            "\n",
            "\n",
            "\n",
            "# --- GRAPH d",
            "\n",
            "set title 'pressure'",
            "\n",
            "set xrange [0:1]",
            "set yrange [0:1]",
            "set size ratio 1",
            "\n",
            "set xtics axis",
            "set xtics format '%.1f'",
            "set xlabel ('x [m]')",
            "\n",
            "set ytics axis",
            "set ytics format '%.1f'",
            "set ylabel ('y [m]')",
            "\n",
            "set cblabel ('p [Pa]')",
            "\n",
            "\n",
            "\n",
            plotp,
            "\n",
            "\n",
            "unset multiplot",
            "pause 50"
        };

        ofstream f;
        string name = "plot_" + to_string(i) + ".gp";
        f.open(name);
        if(!f.is_open()) {
            cout << "error " << i << endl;
            f.close();
        } else {

            for(string s : com)
                f << s << endl;
            f.close();

        }


    }





}
