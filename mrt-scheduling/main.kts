
fun xassert(expr: Boolean) {
    if (!expr)
        throw IllegalStateException()
}

fun error(msg: String = "error") {
    throw IllegalStateException(msg)
}

data class ThetaEntry(val node: Int, var Delta: Double) {
    override fun toString() = "($node, $Delta)"
}

data class FwdArc(val node: Int, val delta: Double)

fun main() {
    val inf = 10e5 

    // problem spec
    /*
    val N = 4
    val lower = listOf<Double>(1.0, 10.0, 0.0, 0.0)
    val upper = listOf<Double>(7.5, 15.0, inf, inf)
    val delta = listOf<Double>(7.5, 10.0, 5.0)
    val isPickup = listOf<Boolean>(true, true, false, false)
    val toPickup = listOf<Int>(-1, -1, 0, 1)
    val theta = (0 .. N-1).map {mutableSetOf<ThetaEntry>()}
    val mrt: Map<Pair<Int, Int>, Double> = mapOf(
            Pair(0, 2) to 18.0,
            Pair(1, 3) to 20.0
    )
    */

    val N = 6
    val lower = listOf<Double>(0.0, 15.0, 0.0, 30.0, 0.0, 0.0)
    val upper = listOf<Double>(inf, inf, inf, inf, inf, inf)
    val delta = listOf<Double>(10.0, 5.0, 5.0, 5.0, 5.0)
    val isPickup = listOf<Boolean>(true, true, false, true, false, false)
    val toPickup = listOf<Int>(-1, -1, 0, -1, 1, 3)
    val theta = (0 .. N-1).map {mutableSetOf<ThetaEntry>()}
    val mrt: Map<Pair<Int, Int>, Double> = mapOf(
            Pair(0, 2) to 15.0,
            Pair(1, 4) to 18.0,
            Pair(3, 5) to 10.0
    )

    // forward pass: tstar[], theta[]
    val tstar = DoubleArray(N)
    val p = mutableListOf<Int>()

    tstar[0] = lower[0]
    if (isPickup[0]) {
        theta[0] += setOf(ThetaEntry(0, 0.0))
        p.add(0)
    }
    for (i in 1 .. N-1) {
        // tstar[]
        if (isPickup[i]) {
            tstar[i] = Math.max(tstar[i-1] + delta[i-1], lower[i])
        } else {
            tstar[i] = tstar[i-1] + delta[i-1]
        }

        // theta[]
        theta[i] += theta[i-1].map { entry ->
            ThetaEntry(entry.node, entry.Delta + delta[i-1])
        }

        if (isPickup[i]) {
            theta[i] += mutableSetOf(ThetaEntry(i, 0.0))
            p.add(i)
        } else
            theta[i].remove(theta[i].first{it.node == toPickup[i]})
    }

    // forward pass: fwdArc[]
    var lastP = -1
    var lastPDelta = 0.0
    val fwdArc: MutableList<FwdArc?> = (0 .. N-1).map {null}.toMutableList()
    for (i in 1 .. N-1) {
        if (isPickup[i]) {
            lastP = i
            lastPDelta = 0.0
        } else if (lastP != -1) {
            lastPDelta += delta[i-1]
            if (lastP != toPickup[i]) {
                mrt[Pair(toPickup[i], i)]?.let {
                    fwdArc[toPickup[i]] = FwdArc(lastP, it - lastPDelta)
                } ?: error("no mrt defined for pair (P=${toPickup[i]}, D=$i)")
            }
        }
    }

    for (i in 0 .. N-1)
        println("t*[$i] = ${tstar[i]} theta = ${theta[i]}")
    println("")

    // backward pass: SP[i]
    val n = p.size
    val SP = DoubleArray(N)

    SP[p[n-1]] = -tstar[p[n-1]]
    for (k in 2 .. n) {
        SP[p[n-k]] = Math.min(SP[p[n-k]], -tstar[p[n-k]])
        for ((pi, piDelta) in theta[p[n-k]]) {
            fwdArc[pi]?.let {
                SP[p[n-k]] = Math.min(SP[p[n-k]], piDelta + it.delta + SP[it.node])
            } ?: error("no fwd arc defined from P=$pi")
        }
    }
     
    // final pass: t[]
    val t = DoubleArray(N)
    for (i in 0 .. n-1)
        println("SP[p$i] = ${SP[p[i]]}")

    t[0] = if (isPickup[0]) -SP[0] else -tstar[0]
    for (i in 1 .. N-1) {
        if (isPickup[i])
            t[i] = -SP[i]
        else
            t[i] = t[i-1] + delta[i-1]
    }

    for (i in 0 .. N-1) {
        print("t[$i] = ${t[i]}")
        if (!isPickup[i])
            print("\tRT = ${t[i]-t[toPickup[i]]} MRT = ${mrt[Pair(toPickup[i], i)]}")
        println("")
    }

}

main()
