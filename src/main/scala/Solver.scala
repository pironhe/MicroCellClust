import scala.collection.mutable.{ListBuffer, PriorityQueue, Map}

import Objective.{buildExprMap, getCoExprMarks, getExpSums, getMarkSum, getMarkers, getMarkersFromPrev, getSampleObj, getObjNoPenalty, getNumberInterstingCells, getNumberNotInterstingCellsInClust}

/*
 * Author: Alexander Gerniers
 */
object Solver {

    /**
      * Find the bicluster maximising the objective using a heuristic search
      * @param m an expression matrix with samples (cells) on the rows and markers (genes) on the columns
      * @param nNeg the maximum percentage of -1 allowed inside the cluster
      * @param kappa a weighting constant for the out-of-cluster expression
      * @param nHeuristic the initial number of top-solutions to consider for expansion at the next level
      * @param maxNbSam the maximum number of samples in the bicluster
      * @param stopNoImprove stop search when no better solution is found after x levels
      * @param minCoExpMark a minimum number of markers that must be expressed in every row of the bicluster
      * @param excl a list of samples that need to be excluded from the search
      * @param verbose enable/disable printing
      * @return - an assignment of samples
      *         - an assignment of markers
      *         - the corresponding objective value
      */
    def findCluster(m: Array[Array[Double]], nP: Array[Array[Int]]= null, nNeg: Double = 0.1, kappa: Double = 1,
                    maxCelWanted: Int= Int.MaxValue, minCelWanted: Int = 0,
                    maxMarkersWanted: Int= Int.MaxValue, minMarkersWanted: Int = 0,
                    nHeuristic: Int = 20,
                    maxNbSam: Int = Int.MaxValue, stopNoImprove: Int = 25, minCoExpMark: Int = 0,
                    excl: List[Int] = List(), verbose: Boolean = true): (List[Int], List[Int], Double, List[List[Int]],List[List[Int]]) = {
        val expr = buildExprMap(m)
        val markSum = getMarkSum(m)
        val nSam = m.length
        val nSamMinusExcl = nSam - excl.length
        var currentKappa = kappa
        var currentNNeg = nNeg
        if (verbose) {
            println("Current search level: 2")
            println("\t " + ((nSamMinusExcl * nSamMinusExcl - nSamMinusExcl) / 2) + " pairs to evaluate")
        }
        val t0 = System.currentTimeMillis
  
  
        val nBestQueue = PriorityQueue[(List[Int], List[Int], Double, Double)]()(Ordering[Double].on(x => -x._3))
        if(nP==null){
            // Evaluation of all the pairs of cells
            
            for (p <- (0 until nSam).combinations(2)) {
                if (excl.length == 0 || !excl.exists(p contains _)) {
                    val (markers, obj, objNoPenalty) = getMarkers(m, p.toList, expr, markSum, nNeg, kappa)
                    nBestQueue += ((p.toList, markers, obj, objNoPenalty))
                    if (nBestQueue.size > nHeuristic) nBestQueue.dequeue
                }
            }
        }else{
            for (p <- nP){
              val(markers,obj,objNoPenalty) = getMarkers(m, p.toList, expr, markSum, nNeg, kappa)
              nBestQueue += ((p.toList, markers, obj, objNoPenalty))
            }
        }
        
        

        val nBest = nBestQueue.dequeueAll.toList.reverse
        //Keep them for later
        val nFirstPairs = nBest
        val nPairs = new ListBuffer[List[Int]]()
        for(t <- nFirstPairs){
          nPairs += t._1
        }
        
        var best = nBest(0)
        var prevLvlNBest = nBest.map(x =>  (x._1, getCoExprMarks(x._1, expr), x._4))
        val t1 = System.currentTimeMillis

        if (verbose) {
            //println("\t cellsPair: " + nPairs.toList)
            println("\t kappa value: "+ currentKappa)
            println("\t mu value: " + currentNNeg)
            println("\t NEW BEST: obj. val.: " + best._3)
            println("\t samples (idx. from 1): " + best._1.map(_ + 1).mkString(", "))
            println("\t nb. markers: " + best._2.length)
            //println("\t obj. val. no penalty: " + best._4)
            println("\t Computation time [s]: " + ((t1 - t0).toDouble / 1000))
        }
        
        val nBestQueueOfLevel = Map[Int,PriorityQueue[(List[Int], List[Int], Double, List[Int], Double)]]()
        val bestSolutionOfLevel = Map[Int,(List[Int], List[Int], Double, List[Int], Double)]()
        val bestSolutionsNbCell = Map[Int, (List[Int], List[Int], Double, Double)]()
        val solutionsSize = ListBuffer[List[Int]]()
        var lvl = 3
        val maxLvl = maxNbSam min nSam
        var noImprove = 0
        var retVal = (List[Int](), List[Int](), 0.0, List[List[Int]]())
        var finished = false
        var restartComputing = false
        var newSolFound = false
        var numberGeneOk = false
        var numberCellOk = false
        var numberComputation = 0
        var modifNNeg = 0
        while ( !finished && (!numberGeneOk || !numberCellOk) && lvl <= maxLvl) {
            
            if(restartComputing){ 
                println("NUMBER COMPUTATION: " + numberComputation)
                //Compute markers, obj val and objval no neg again with new kappa and/or mu
                //Keep the pairs of cell selected
                //recup pair
                //recalculer markers
                //ajouter à nBestQueue
                //quand recompute reinit lvl, noimprove
                if(verbose){
                    println("\t value of mu : " + currentNNeg)
                    println("\t value of kappa : " + currentKappa)
                }
                val nBestQueue = PriorityQueue[(List[Int], List[Int], Double, Double)]()(Ordering[Double].on(x => -x._3))
                for(p <- nFirstPairs.indices){
                    val samples = nFirstPairs(p)._1
                    val (markers, obj, objNoPenalty) = getMarkers(m, samples, expr, markSum, currentNNeg, currentKappa)
                    nBestQueue += ((samples, markers, obj, objNoPenalty))
                    if (nBestQueue.size > nHeuristic) nBestQueue.dequeue
                }
                //recomputation done
                val nBest = nBestQueue.dequeueAll.toList.reverse
                best = nBest(0)
                prevLvlNBest = nBest.map(x =>  (x._1, getCoExprMarks(x._1, expr), x._4))
                noImprove = 0
                println("Restarted Computing")
                if (verbose) {
                    println("\t NEW BEST: obj. val.: " + best._3)
                    println("\t samples (idx. from 1): " + best._1.map(_ + 1).mkString(", "))
                    println("\t nb. markers: " + best._2.length)
                    println("\t obj. val. no penalty: " + best._4)
                }
                restartComputing = false
            }
            if (verbose) {
                println("Current search level: " + lvl)
            }
            val t0 = System.currentTimeMillis
            val nBestQueue = PriorityQueue[(List[Int], List[Int], Double, List[Int], Double)]()(Ordering[Double].on(x => -x._3))

            for (p <- prevLvlNBest.indices) {
                val prevExpSums = getExpSums(m, prevLvlNBest(p)._1, nNeg)
                var toExclude = excl ++ prevLvlNBest(p)._1
                for (pp <- 0 until p) {
                    val setDif = prevLvlNBest(pp)._1.filterNot(prevLvlNBest(p)._1.toSet)
                    if (setDif.length == 1) toExclude ++= setDif
                }
                for (i <- 0 until nSam if !toExclude.contains(i)) {
                    val samples = (i :: prevLvlNBest(p)._1).sorted
                    val coExp = (prevLvlNBest(p)._2.toSet intersect expr(i).toSet).toList
                    if (coExp.length >= minCoExpMark) {
                        val (markers, obj) = getMarkersFromPrev(m, samples, i, markSum, prevExpSums, currentNNeg, currentKappa)
                        val objNoPenalty = getObjNoPenalty(m , samples, markers)
                        nBestQueue += ((samples, markers, obj, coExp, objNoPenalty))
                        if (nBestQueue.size > nHeuristic) nBestQueue.dequeue
                    }
                }
            }

            if (nBestQueue.size > 0) {
                nBestQueueOfLevel += (lvl -> nBestQueue)
                val nBest = nBestQueue.dequeueAll.toList.reverse
                prevLvlNBest = nBest.map(x => (x._1, x._4, x._5))
                var bestOfLvlOption = bestSolutionOfLevel.get(lvl)
                if (bestOfLvlOption!=None){
                    var bestOfLvl = bestOfLvlOption.get
                    if(nBest(0)._3 > bestOfLvl._3 && nBest(0)._2.length <= maxMarkersWanted && nBest(0)._2.length >= minMarkersWanted){
                        bestSolutionOfLevel += (lvl -> nBest(0))
                    }
                  
                }
                if (nBest(0)._3 > best._3) {
                    newSolFound = true
                    best = (nBest(0)._1, nBest(0)._2, nBest(0)._3, nBest(0)._5)
                    bestSolutionsNbCell += (best._1.length -> best)
                    noImprove = 0
                    if (verbose) {
                        println("\t NEW BEST: obj. val.: " + best._3)
                        println("\t samples (idx. from 1): " + best._1.map(_ + 1).mkString(", "))
                        println("\t nb. markers: " + best._2.length)
                        //println("\t obj. val. no penalty: " + best._4)
                    }
                } else {
                    noImprove += 1
                    if (noImprove >= stopNoImprove) {
                        /*var cellMeanObjValue= best._4 / best._1.length
                        var nbInterestingCells = getNumberInterstingCells(m, best._1, best._2, cellMeanObjValue)
                        println("\t INSIDE cellMeanObjValue: "+ cellMeanObjValue)
                        println("\t nbInterestingCells: "+ nbInterestingCells)*/
                        finished = true
                        if (verbose) {
                            println("\t No improvement after " + stopNoImprove + " levels: search stopped")
                        }
                    }
                }
                lvl += 1
            } else {
              finished = true
            }
            if(finished){
                numberComputation +=1
                if(numberComputation==1){
                    retVal = (best._1,best._2,best._3, nPairs.toList)
                    println("\t Inside NEW BEST Solution, obj val: " + best._3)
                }else{
                    if(retVal._3 < best._3 && (best._2.length >= 10 && best._2.length <=30) 
                    || ((retVal._2.length < 10 || retVal._2.length > 30) && 
                    (best._2.length >= 10 && best._2.length <=30))){
                        // If the new solution has a better obj value and has a good cluser size
                        // or it has a good cluster size while the current one hasn't
                        retVal = (best._1,best._2,best._3, nPairs.toList)
                        println("\t Inside NEW BEST Solution, obj val: " + best._3)
                    }
                }
                
                println("\t RETVAL currnet: " + retVal)
                solutionsSize += List(best._1.length, best._2.length, best._3.toInt)
                if(numberComputation==10){
                    finished = true
                }else{
                    finished= false
                    restartComputing = true
                    lvl = 3
                    if(best._2.length > 30 && best._1.length > 5){ //need to increase kappa
                        currentKappa = 1.2*currentKappa
                        modifNNeg =0
                    }else if(best._2.length > 30 && best._1.length <= 5){//need to decrease kappa solution too small
                        currentKappa = 0.8*currentKappa
                        modifNNeg =0
                    }else if (best._2.length < 10){ //need to decrease kappa
                        currentKappa = 0.8*currentKappa
                        modifNNeg =0
                    }else{ // TUNING OF MU HERE
                        var cellMeanObjValue= best._4 / best._1.length
                        var nbInterestingCells = getNumberInterstingCells(m, best._1, best._2, cellMeanObjValue)
                        println("\t INSIDE cellMeanObjValue: "+ cellMeanObjValue)
                        println("\t nbInterestingCells: "+ nbInterestingCells)
                        if(nbInterestingCells > 0 && modifNNeg==0){
                            var percentageModif = 1.0 + (nbInterestingCells.toDouble / best._1.length.toDouble)
                            currentNNeg = currentNNeg * percentageModif
                            modifNNeg +=1
                            println("\t new nneg value : " + currentNNeg)
                            println("\t percentage modif: " + percentageModif)
                        }else{
                          var nbNotInterestingCells = getNumberNotInterstingCellsInClust(m, best._1, best._2, cellMeanObjValue)
                          println("\t nbNotInterestingCells: "+ nbNotInterestingCells)
                          if(nbNotInterestingCells > 0 && modifNNeg < 2 ){
                            modifNNeg += 1
                            var percentageModif = 1.0 - ( nbNotInterestingCells.toDouble / best._1.length.toDouble)
                            currentNNeg = currentNNeg * percentageModif
                            println("\t percentage modif: " + percentageModif)
                          }else{
                            finished = true
                          }
                        }
                        
                    }
                }
            }
            
            /*if(lvl>maxCelWanted){ // No need to go that far in the search
                // Check number of cells in best solution to choose if need to increase or decrease nNeg
                // best = (samples, markers, obj, objNoPenalty)
                if( best._1.length < minCelWanted ){ //Need to increase nNeg
                    var cellMeanObjValue= best._4 / best._1.length
                    var nbInterestingCells = getNumberInterstingCells(m, best._1, best._2, cellMeanObjValue)
                    currentNNeg *= 1.1
                    restartComputing = true
                    println("\t INSIDE cellMeanObjValue: "+ cellMeanObjValue)
                    println("\t nbInterestingCells: "+ nbInterestingCells)
                    lvl = 3
                }else if( best._1.length > maxCelWanted ){ //Need to decrease nNeg
                    lvl = 3
                    currentNNeg = currentNNeg*0.9
                    restartComputing = true
                }else{
                    numberCellOk = true
                    println("\t Number Cells ok")
                }
            }
            if(finished){
                if( best._2.length < minMarkersWanted ){ //Need to reduce kappa
                    currentKappa *= 0.9
                    restartComputing = true
                    lvl = 3
                    finished = false
                }else if( best._2.length > maxMarkersWanted ){ //Need to increase kappa
                    currentKappa *= 1.1
                    restartComputing = true
                    lvl = 3
                    finished = false
                }else{
                    numberGeneOk = true
                    println("\t Number Markes ok")
                }
            }*/
            
            
            
            val t1 = System.currentTimeMillis
            if (verbose) {
                println("\t Computation time [s]: " + ((t1 - t0).toDouble / 1000))
            }
        }
        
        //Return the best solution which has a number of cell and markers inside the wanted area
        
        return (retVal._1,retVal._2,retVal._3,retVal._4,solutionsSize.toList)
    }

    /**
      * Get the objective value of each pair of samples
      * @param m an expression matrix with samples (cells) on the rows and markers (genes) on the columns
      * @param kappa a weighting constant for the out-of-cluster expression
      * @return the objective values in decreasing order
      */
    def evaluatePairs(m: Array[Array[Double]], kappa: Double = 1): List[(List[Int], Double)] = {
        val expr = buildExprMap(m)
        val markSum = getMarkSum(m)
        val nSam = m.length

        val lb = ListBuffer[(List[Int], Double)]()

        for (p <- (0 until nSam).combinations(2)) {
            val (_, obj, objNoPenalty) = getMarkers(m, p.toList, expr, markSum, 0, kappa)
            lb += ((p.toList, obj))
        }

        return lb.toList.sortBy(-_._2)
    }

}