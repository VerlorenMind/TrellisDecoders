gnuplot -persist <<-EOFMarker
	set terminal eps enhanced font "Times-New-Roman,10" linewidth 3
	set output "$1_FER_.eps"
	set datafile separator ","
	set logscale y 10
	set grid ytics lc rgb "#bbbbbb" lw 1 lt 0
	set grid xtics lc rgb "#bbbbbb" lw 1 lt 0
	set grid mytics lc rgb "#bbbbbb" lw 1 lt 0
	set ylabel "{/*1.5 Word error probability}"
	set xlabel "{/*1.5 E_b/N_0, dB}"
        set format y "10^{%L}"
	plot "$1" using 1:2 with lines ti "BEAST", "$1" using 1:3 with lines ti "BSD"
        unset logscale y
	unset yrange
        unset format y
	set output "$1_sum_.eps"
	set ylabel "{/*1.5  Summation operations}"
	set xlabel "{/*1.5 E_b/N_0, dB}"
        plot "$1" using 1:5 with lines ti "BEAST", "$1" using 1:7 with lines ti "BSD"
        set output "$1_comp_.eps"
	set ylabel "{/*1.5 Comparison operations}"
	set xlabel "{/*1.5 E_b/N_0, dB}"
        plot "$1" using 1:6 with lines ti "BEAST", "$1" using 1:8 with lines ti "BSD"
EOFMarker
