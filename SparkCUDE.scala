package com.island.SparkTest

import java.util
import javax.naming.ConfigurationException

import cec2010._
import cer2003.{Common, FitnessFunctions}
import com.island.SparkStrategies
import com.island.SparkStrategies._
import org.apache.spark.{SparkContext, SparkConf}
import org.apache.spark.SparkContext._

import java.util.logging.Level;
import java.util.logging.Logger;
//������
import java.util.{HashMap => JHashMap}
import scala.collection.mutable.HashMap
import scala.collection.immutable.HashMap

object SparkCUDE {  
  var path = "/resultPath"
  def main(args: Array[String]) {      
      val nRepeatTest:Int = 25;// 
      val dimension:Int = 1000;//
      val rounds:Int = 30;//
      val islandCount:Int = 5;//
      val generationsPerRound:Int = 1000;//
    
      val defaultPopulationSize:Int = 100/islandCount;
      val migratingIndividuals:Int = 15;
      //����ÿ�ζ���ִ�е�����ֵ
      val win = new Array[Double](nRepeatTest);

      //��ֹʱ��
      var StartTime: Long = 0
      var EndTime: Long = 0
      val T = new Array[Long](nRepeatTest);

      //���ƿ���
      val deleteHistory = true;
      val verbose = false;
      val writeHistory = true;
      val runFirstTest = true;
      //-------------variables-------------      
      var deaf: txjDEDeaf=null
      for(f<-1 to 20){//20
        //------------CER2003--- Test_APP-------------        
        for(r <- 0 to nRepeatTest-1){          
          f match {
            case 1=>              
              val fitFunction=new F1()               
              val derand1bin= new CUDE(fitFunction,dimension,defaultPopulationSize)
              deaf = new txjDEDeaf(path,fitFunction, dimension,derand1bin , generationsPerRound, islandCount, defaultPopulationSize, migratingIndividuals)
            case 2=>              
              val fitFunction=new F2()
              val derand1bin= new CUDE(fitFunction,dimension,defaultPopulationSize)
              deaf = new txjDEDeaf(path,fitFunction, dimension,derand1bin , generationsPerRound, islandCount, defaultPopulationSize, migratingIndividuals)
            case default=>deaf=null;System.out.println("no fitnessfunctions\n");System.exit(0);
          }         
          val Fs=new Array[Double](islandCount)
          val CRs=new Array[Double](islandCount)
          for (t <- 0 to islandCount-1) { //algorithm parameters
            Fs.update(t,0.5)
            CRs.update(t,0.9)
          }
          //ÿ�������ղ�ͬ��F��CR���Ż�����
          for (i <- 0 to islandCount-1) { //set parameters to algorithms            
            val subpopParameters = new JHashMap[String, Double]();            
            subpopParameters.put("F", Fs.apply(i))
            subpopParameters.put("CR", CRs.apply(i))            
            deaf.addSubpopulationsConfig(i,classOf[DERand1bin], subpopParameters);
          }
          deaf.setVerbose(verbose);          
          //���������˽��
          deaf.setTopology(classOf[Topologies.Ring_1_2]);
          //�����滻����
          deaf.setImmigrationMethod(SiPDEPopulation.acceptImmigrantsMethod.REPLACE_WORST);
          deaf.setEmigrationMethod(SiPDEPopulation.expelEmigrantsMethod.EXPEL_BEST);
          try {
            System.out.println("********************* Test " + r + " *********************")
            StartTime=System.currentTimeMillis()
            //��ʼ����Ⱥ
            var rdd= deaf.createRandomPopulation(SparkDE.sc)
            //��Ⱥ����
            var winner:(Int, SiPDEPopulation) = deaf.run(rdd,rounds);            
            EndTime=System.currentTimeMillis()            
            T.update(r,EndTime - StartTime)            
            System.out.println("����ʱ��=="+T.apply(r))
            val tempvalue="����ʱ��=="+T.apply(r)+"And the winner is : " +winner._1+"And the winner is : " +winner._2.getBestIndividual.getFitness+"\n"
            //��tempvalue�����path��
	    Common.appendMethodB(path, tempvalue)            
          } catch {
            case ex:ConfigurationException =>{
              Logger.getLogger(SiPDETest.getClass.getName()).log(Level.SEVERE, null, ex);
            }            
          }
        }
      }
    }
}
object SparkDE {
  val sparkConf = new SparkConf()
			.setAppName("SparkDECC")
			.setMaster("local")
			.set("spark.executor.memory","2g")  
  val sc = new SparkContext(sparkConf)  
}
