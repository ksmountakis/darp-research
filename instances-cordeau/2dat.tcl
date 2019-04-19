#!/usr/bin/tclsh

proc prob_load {fd} {
	set prob [dict create]

	# read first line
	# <num> <num> <max veh ride minutes Tk> <veh capacity Qk> <max user ride minutes Li>
	if {[gets stdin line] < 0} {
		error "No input"
	}
	dict set prob m [lindex $line 0]
	dict set prob T [lindex $line 2]
	dict set prob Q [lindex $line 3]
	dict set prob L [lindex $line 4]

	# read remaining lines
	# <i> <x coord xi> <y coord yi> <service time at node di> <num passengers qi> <earliest time ei> <latest time li>
	# note: traveling time from i to j = eucledian distance based on (xi, yi) and (xj, yj)

	dict set prob n 0
	while {[gets $fd line] >= 0} {
		lassign $line i xi yi di qi ei li
		dict set prob x$i $xi
		dict set prob y$i $yi
		dict set prob d$i $di
		dict set prob q$i $qi
		dict set prob e$i $ei
		dict set prob l$i $li

		for {set j 0} {$j <= [dict get $prob n]} {incr j} {
			set xj [dict get $prob x$j]
			set yj [dict get $prob y$j]
			set dx [expr pow($xj - $xi,2)]
			set dy [expr pow($yj - $yi,2)]
			dict set prob t$i:$j [expr sqrt($dx+$dy)]
			dict set prob t$j:$i [expr sqrt($dx+$dy)]
		}

		dict incr prob n
	}

	return $prob
}

proc ampl_matrix {key dictionary dimX dimY} {
	puts -nonewline "param $key :\n  "
	for {set i 0} {$i < $dimX} {incr i} {
		puts -nonewline [format "%3d " [expr $i+1]]
	}
	puts ":="

	for {set i 0} {$i < $dimY} {incr i} {
		puts -nonewline "[expr $i+1] "
		for {set j 0} {$j < $dimX} {incr j} {
			set tij [dict get $dictionary ${key}${i}:${j}]
			puts -nonewline [format "%-2.3f " $tij]
		}
		if {$i == [expr $dimX-1]} {
			puts -nonewline ";"
		}
		puts ""
	}
}

proc ampl_scalar {key value} {
	puts "param $key := $value;"
}

proc ampl_vector {key dictionary dimX} {
	puts "param $key := "
	for {set i 0} {$i < $dimX} {incr i} {
		puts -nonewline "[expr $i+1] [dict get $dictionary ${key}$i]"
		if {$i == [expr $dimX-1]} {
			puts -nonewline ";"
		}
		puts ""
	}
}

proc prob_ampl {p} {
	set n [dict get $p n]
	ampl_scalar "n" [expr ($n-2)/2]
	ampl_scalar "m" [dict get $p m]
	ampl_scalar "L" [dict get $p L]
	ampl_scalar "Q" [dict get $p Q]
	ampl_scalar "T" [dict get $p T]
	ampl_matrix "t" $p $n $n
	ampl_vector "d" $p $n
	ampl_vector "q" $p $n
	ampl_vector "e" $p $n
	ampl_vector "l" $p $n
}

set p [prob_load stdin]
puts [prob_ampl $p]
