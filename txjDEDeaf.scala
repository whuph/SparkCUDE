package com.island.SparkStrategies

import java.{util, lang}

import com.island.SparkStrategies.Topologies.Ring
import com.util.{SiPDEPartitioner, IteblogPartitioner, newInstanceObject}

import org.apache.spark.rdd.RDD

import org.apache.spark._
import org.apache.spark.SparkContext

import org.apache.spark.SparkContext._
import org.apache.spark.{SparkConf,SparkContext}
import java.lang.reflect.Constructor
//重命名
import java.util.{HashMap => JHashMap}
import scala.collection.mutable.{ArrayBuffer, HashMap}
import scala.collection.immutable.HashMap
import cec2010.Function
;
/**
 *  16-11-12.
 */
class txjDEDeaf extends java.io.Serializable{
    private var nMapTasks: Int = 3
    private[SparkStrategies] var VERBOSE: Boolean = false
    private[SparkStrategies] var WRITEHISTORY: Boolean = false
    private var defaultPopSize: Int = 0
    private[SparkStrategies] var islands: Int = 0
    private[SparkStrategies] var migrationRounds: Int = 0
    private[SparkStrategies] var round: Integer = 0
    var dimensions: Int = 0
    var migratingPopSize: Int = 1
    private[SparkStrategies] var generations: Int = 0
    private[SparkStrategies] var subPopulationSizes: util.ArrayList[Int] =new util.ArrayList[Int]{}
    private[SparkStrategies] var defaultPopulationAlgorithm: Class[_] = null
    private[SparkStrategies] var subPopulationAlgorithms: util.ArrayList[Class[_]] = new util.ArrayList[Class[_]]
    private var function: Class[_] = null
    private var algorithm:Algorithm = null;
    private var fInstance: Function = null
    private var configuredSubpopulations: Int = 0
    private var topology: Class[_] = classOf[Topology]
    private var immigrationMethod: SiPDEPopulation.acceptImmigrantsMethod = null
    private var emigrationMethod: SiPDEPopulation.expelEmigrantsMethod = null
    private[SparkStrategies] var DELETE_HISTORY: Boolean = true
    private[SparkStrategies] var RANDOM_ALGS_AND_PARAMS: Boolean = false
    private[SparkStrategies] var NO_OF_ALGORITHMS: Int = 5
    //起止时间
    var StartTime: Long = 0
    var EndTime: Long = 0
    var path=""    
    private var configDouble:util.ArrayList[JHashMap[String,Double]]=new util.ArrayList[JHashMap[String, Double]]()
    private val configClass: JHashMap[String, String] = new JHashMap[String, String]()    
    private val genValue:JHashMap[Int, Double] = new JHashMap[Int, Double]()
    var strArrayVar=ArrayBuffer((0,0.0))    
    def this(temppath:String,fitFunction:  Function, D: Int, algorithm:  Algorithm, generationsPerRound: Int, islandCount: Int, populationSize: Int, migratingIndividuals: Int) {
        this()
        this.path=temppath
        generations = generationsPerRound
        islands = islandCount
        defaultPopSize = populationSize
        nMapTasks = islandCount
        migratingPopSize = migratingIndividuals
        this.algorithm=algorithm
        dimensions = D        
        fInstance=fitFunction
        round = 0
    }
    private def configureInitialization(sc: SparkContext):RDD[(Int, SiPDEPopulation)] ={
        val rdd=sc.parallelize(0 to islands-1, islands).map(i=>{
            val popSize: Int = defaultPopSize  //每个岛中子群的数目            
            val population: SiPDEPopulation = new SiPDEPopulation
            population.clear
            population.setFitnessFunction(fInstance)//设置函数；
            population.addRandomIndividuals(popSize, dimensions)
            population.resetIndividualsPopulation
            population.setKey(i)//设置岛的编号
            (i,population)
        }).cache()
        rdd
    }
    //从所有岛中寻找最优的个体
    private def findBestIndividual(rdd:RDD[(Int,SiPDEPopulation)]): (Int, SiPDEPopulation) = {
        val minrdd= rdd.reduce((x,y)=>{
            if(x._2.getBestIndividual.getFitness<y._2.getBestIndividual.getFitness){
                (x._1,x._2)
            }else{
                (y._1,y._2)
            }
        })        
        minrdd
    }

    private def printBestIndividual(winner: (Int, SiPDEPopulation)) {
        if (winner == null) {
            System.out.println("---- Unable to print best individual !")
        }
        else {
            System.out.println("---- Best so far : \n" + winner._2.getBestIndividual.getFitness)
        }
    }

    def createRandomPopulation(sc: SparkContext): RDD[(Int, SiPDEPopulation)] = {
        val rdd=configureInitialization(sc)        
        val bestInd: (Int, SiPDEPopulation) = findBestIndividual(rdd)        
        rdd
    }

    var bestInd: (Int, SiPDEPopulation) = null
    def run(rdd:RDD[(Int, SiPDEPopulation)],round:Int):((Int, SiPDEPopulation))={
        //最优个体        
        var generationsPerRound: Int = 0
        var oldRdd=rdd 
        //迁移的拓扑结构
        var  topology:Topology =new Ring(islands);       
        var tempvalue=""
        val tempvalue1 = " F " + (fInstance.getClass.getName) +"  run " +  "  " + " Slices " + "\n"
        Common.appendMethodB(path, tempvalue1)
        //var genValue=null
        for (i <-0 to round-1) {
            //迁入的种群
            val rdd1= oldRdd.map(r=>{
                algorithm.setPopulation(r._2)
                r._2.resetIndividualsPopulation()
                algorithm.setParameter("F",configDouble.get(r._1).get("F"))
                algorithm.setParameter("CR",configDouble.get(r._1).get("CR"))
                generationsPerRound=generations//每个子种群中的进化代数
                for(j<-0 to generationsPerRound-1){
                    val sipdeindividuals:SiPDEIndividuals= algorithm.generation();
                }                
                val pop=algorithm.getPopulation
                (r._1,pop)
            })            
            //迁出的种群  提出最优个体
            val emigrants:RDD[(Int, SiPDEPopulation)]=rdd1.map(r=>{
                val emig=r._2.getEmigrants(migratingPopSize,emigrationMethod)
                emig.resetIndividualsPopulation();
                (r._1,emig)
            })
            //替换最差个体
            val partitioner=new SiPDEPartitioner(emigrants,topology,islands,immigrationMethod)
            val rdd2=rdd1.map(r=>{
                (r._2,r._1)
            }).partitionBy(partitioner).cache()
            //将迁出的种群传递到分区函数中
            oldRdd=rdd2.map(r=>{
                (r._2,r._1)
            })
            StartTime=System.currentTimeMillis()
            var tempbestInd = findBestIndividual(oldRdd)
            tempvalue += (i+1)+" "+tempbestInd._2.getBestIndividual.getFitness
            EndTime=System.currentTimeMillis()
            tempvalue+= " "+(EndTime - StartTime)+" \n"
        }
        //输出最优的值
        bestInd = findBestIndividual(oldRdd)        
        Common.appendMethodB(path, tempvalue)
        printBestIndividual(bestInd)
        bestInd
    }

    def getRound: Int = {
        round
    }

    def setRound(nr: Int) {
        round = nr
    }

    def getGenValue(): JHashMap[Int, Double] ={
        genValue
    }

    /**
     * Sets algorithms according to their ratings...
     */
    def setRatedAlgorithms {
    }
    def setTopology(top: Class[_ <: Topology]) {
        topology = top
    }    
    def setImmigrationMethod(met: SiPDEPopulation.acceptImmigrantsMethod) {
        immigrationMethod = met
    }
    def setEmigrationMethod(met: SiPDEPopulation.expelEmigrantsMethod) {
        emigrationMethod = met
    }
    def setVerbose(verb: Boolean) {
        VERBOSE = verb
    }
    
    def addSubpopulationsConfig(index:Int,alg:Class[_ <: Algorithm], subpopParameters:JHashMap[String, Double]) {
        addSubpopulationsConfig(index,alg, defaultPopSize, subpopParameters);
    }
    //addSubpopulationsConfig
    def addSubpopulationsConfig(index:Int,alg:Class[_ <: Algorithm], popSize:Int, subpopParameters:JHashMap[String, Double]){
        //将参数F、CR存放在相应的岛中
        configDouble.add(index,subpopParameters)        
        //存放优化器
        configClass.put("algorithm_"+index,alg.toString)
        //种群大小
        configClass.put("popSize_"+index,popSize.toString)

    }
}