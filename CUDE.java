package com.island.SparkStrategies;

import cec2010.Function;

import com.java.tools.RandSet;
import com.util.rnd;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

public class CUDE extends Algorithm {
	static protected double F = 0.1;
	static protected double CR = 0.1;

	static protected double lowF = 0.5;
	static protected double uppF = 0.8;
	static protected double lowCR = 0.1;
	static protected double uppCR = 0.9;
	static protected double std = 0.1;
	static protected int sn = 6;
	protected int NP = popSize;
	protected double[][] dictsucc = new double[NP][sn];
	protected double[][] dicttota = new double[NP][sn];
	int Q=6;
	int J=6;
	int[][] temp_layout=new int[Q+1][J];
	int[][] uniform_layout=new int[Q][J];

	public CUDE() {

	}

	public CUDE(Function f, int D_, int popSize_) {
		dimensions = D_;
		popSize = popSize_;
		function = f;
		population.setFitnessFunction(f);
		minPopulationSize = 5;
		lowF = 0.5;
		uppF = 0.8;
		lowCR = 0.1;
		uppCR = 0.9;
		std = 0.1;
		Q=6;
		J=6;
		NP = popSize;
		sn = 6;
		dictsucc = new double[NP][sn];
		dicttota = new double[NP][sn];
		ones(dictsucc);
		ones(dicttota);
		temp_layout=Uniform_array_generate(Q+1, J);
		//uniform_layout=;
		getSubUnifor(uniform_layout,temp_layout,Q);
	}

	@Override
	public SiPDEIndividuals generation() throws IOException {
		// TODO Auto-generated method stub
		try {
			if (popSize < minPopulationSize) {
				throw new Exception("popSize can't be smaller than " + minPopulationSize + "");
			}
		} catch (Exception ex) {
			Logger.getLogger(DERand1bin.class.getName()).log(Level.SEVERE, null, ex);
		}

		SiPDEIndividuals rand1, rand2, rand3;
		RealVector noisy;
		RealVector trial;

		double[] active;
		//System.out.println("popSize="+popSize);
		double trialFitness, activeFitness;
		for (int ind = 0; ind < popSize; ind++) {
			trial = new ArrayRealVector(dimensions);
			noisy = new ArrayRealVector(dimensions);
			active = population.get(ind).getGeno().toArray();
			activeFitness = population.get(ind).getFitness();
			int[] randIndex= RandSet.randomArray1(0, popSize - 1, 3, ind);			
			rand1 = population.get(randIndex[0]);
			rand2 = population.get(randIndex[1]);
			rand3 = population.get(randIndex[2]);

			double[] sr = division(dictsucc[ind], dicttota[ind]);

			double[] P = division(sr, sum(sr));
			int scheme = Roulette(P);
			SiPDEIndividuals best = null;
			double randF=0;
			switch (scheme) {
				case 1:
					F=lowF + Math.random()*std;
					CR=lowCR+Math.random()*std;
					for (int j = 0; j < dimensions; j++) {
						double pom = (rand1.getGene(j) - rand2.getGene(j)) * F + rand3.getGene(j);
						if (pom < function.getMin()|| pom > function.getMax()) {
							pom =function.getMin()+(Math.random() * ((function.getMax() - function.getMin()) + 1)); //function.getRandomValueInDomains(j);
						}
						noisy.addToEntry(j, pom);
						// Create trial Individual from noisy and active
						trial.addToEntry(j, (Math.random() < CR) ? noisy.getEntry(j) : active[j]);
					}
					break;
				case 2:
					F=uppF + Math.random()*std;
					CR=uppCR+Math.random()*std;
					for (int j = 0; j < dimensions; j++) {
						double pom = (rand1.getGene(j) - rand2.getGene(j)) * F + rand3.getGene(j);
						if (pom < function.getMin()|| pom > function.getMax()) {
							pom = function.getMin()+(Math.random() * ((function.getMax() - function.getMin()) + 1));//function.getRandomValueInDomains(j);
						}
						active[j];
						noisy.addToEntry(j, pom);
						trial.addToEntry(j, (Math.random() < CR) ? noisy.getEntry(j) : active[j]);
					}
					break;
				case 3:
					F=lowF + Math.random()*std;
					CR=lowCR+Math.random()*std;
					best = population.getBestIndividual();
					for (int j = 0; j < dimensions; j++) {
						double pom = (rand1.getGene(j) - rand2.getGene(j)) * F + + best.getGene(j);
						if (pom < function.getMin()|| pom > function.getMax()) {
							pom = function.getMin()+(Math.random() * ((function.getMax() - function.getMin()) + 1));
						}
						noisy.addToEntry(j, pom);						
						trial.addToEntry(j, (Math.random() < CR) ? noisy.getEntry(j) : active[j]);
					}
					break;
				case 4:
					F=uppF + Math.random()*std;
					CR=uppCR+Math.random()*std;
					best = population.getBestIndividual();
					for (int j = 0; j < dimensions; j++) {//rand3.getGene(j)
						double pom = (rand1.getGene(j) - rand2.getGene(j)) * F + + best.getGene(j);
						if (pom < function.getMin()|| pom > function.getMax()) {
							pom = function.getMin()+(Math.random() * ((function.getMax() - function.getMin()) + 1));//function.getRandomValueInDomains(j);
						}						
						noisy.addToEntry(j, pom);
						trial.addToEntry(j, (Math.random() < CR) ? noisy.getEntry(j) : active[j]);
					}
					break;
				case 5:
					F=lowF + Math.random()*std;
					CR=lowCR+Math.random()*std;
					randF=Math.random();
					for (int j = 0; j < dimensions; j++) {//rand3.getGene(j)
						double temp =population.get(ind).getGene(j); //(rand1.getGene(j) - rand2.getGene(j)) * F + + best.getGene(j);
						double pom =temp+randF*(rand1.getGene(j) - temp)
								+F*(rand2.getGene(j) - rand3.getGene(j));
						if (pom < function.getMin()|| pom > function.getMax()) {
							pom = function.getMin()+(Math.random() * ((function.getMax() - function.getMin()) + 1));//function.getRandomValueInDomains(j);
						}						
						noisy.addToEntry(j, pom);
						trial.addToEntry(j, (Math.random() < CR) ? noisy.getEntry(j) : active[j]);
					}
					break;
				case 6:
					F=uppF + Math.random()*std;
					CR=uppCR+Math.random()*std;
					randF=Math.random();
					for (int j = 0; j < dimensions; j++) {//rand3.getGene(j)
						double temp =population.get(ind).getGene(j); //(rand1.getGene(j) - rand2.getGene(j)) * F + + best.getGene(j);
						double pom =temp+randF*(rand1.getGene(j) - temp)
								+F*(rand2.getGene(j) - rand3.getGene(j));
						if (pom < function.getMin()|| pom > function.getMax()) {
							pom = function.getMin()+(Math.random() * ((function.getMax() - function.getMin()) + 1));//function..getRandomValueInDomains(j);
						}						
						noisy.addToEntry(j, pom);
						trial.addToEntry(j, (Math.random() < CR) ? noisy.getEntry(j) : active[j]);
					}
					break;
			}
			trialFitness = function.compute(trial.toArray());
			//修改dictsucc、dicttota数组中的值
			if(trialFitness<population.get(ind).getFitness()){
				dictsucc[ind][scheme-1]=dictsucc[ind][scheme-1]+1;
			}
			dicttota[ind][scheme-1]=dicttota[ind][scheme-1]+1;
			if (trialFitness < activeFitness) {
				SiPDEIndividuals indiv = new SiPDEIndividuals(population, trial);
				indiv.setFitness(trialFitness);
				population.set(ind, indiv);
			}

			if (population.get(ind).getFitness() < bestFitness) {
				bestIndividual = population.get(ind);
				bestFitness = bestIndividual.getFitness();
			}
		}		
		//Uniform local search 
		int[] bid=rnd.randomArray(0, NP-1, NP);
		double[][] U = new double[NP][dimensions];		
		List<int[]> subvectors =Uniform_local_search(U, population.get(bid[0]).getGeno().toArray(),
				population.get(bid[1]).getGeno().toArray()
				, dimensions, uniform_layout, Q, J);
		double[] fit_u=new double[Q];
		for(int j=0;j<Q;j++){
			fit_u[j]=function.compute(U[j]);
		}
		double[] minfit=min(fit_u);
		if(bid[0]>minfit[0]){			
			////更改该种群中个体bid(0)
			SiPDEIndividuals indiv = new SiPDEIndividuals(population,new ArrayRealVector(U[(int)minfit[1]]));
			indiv.setFitness(minfit[0]);
			population.set(bid[0], indiv);
		}else if(bid[1]>minfit[0]){
			////更改该种群中个体bid(1)
			SiPDEIndividuals indiv = new SiPDEIndividuals(population, new ArrayRealVector(U[(int)minfit[1]]));
			indiv.setFitness(minfit[0]);
			population.set(bid[1], indiv);
		}
		bestFitness=population.getBestIndividual().getFitness();			
		return bestIndividual;
	}

	@Override
	public SiPDEPopulation getPopulation() {
		// TODO Auto-generated method stub
		return population;
	}

	@Override
	public void newRound() {
		// TODO Auto-generated method stub

	}

	public void setParameter(String configName, double value) {
		// TODO Auto-generated method stub
		if (configName.equals("F")) {
			F = value;
		} else if (configName.equals("CR")) {
			CR = value;
		}
	}

	/**
	 * @param sr
	 * @return sum
	 */
	public double sum(double[] sr) {
		double value = 0;
		for (int i = 0; i < sr.length; i++) {
			value += sr[i];
		}
		return value;
	}

	/**
	 * @param sr
	 * @return 求一个数组的最小值，返回一个数组[最小值，小标]
	 */
	public double[] min(double[] sr) {
		double[] value = new double[2];
		double min=sr[0];
		int k=0;
		for (int i = 1; i < sr.length; i++) {
			//value += sr[i];
			if(sr[i]<min){   // 判断最小值
				min=sr[i];
				k=i;
			}
		}
		value[0]=min;
		value[1]=k;
		return value;
	}

	/**
	 * @param mol
	 * @param Den
	 * @return 两个数组对应值相除
	 */
	public double[] division(double[] mol, double[] Den) {
		double[] value = new double[mol.length];
		for (int i = 0; i < mol.length; i++) {
			value[i] = mol[i] / Den[i];
		}
		return value;
	}

	/**
	 * @param mol
	 * @param Den
	 * @return 两个数组对应值相除
	 */
	public double[] division(double[] mol, double Den) {
		double[] value = new double[mol.length];
		for (int i = 0; i < mol.length; i++) {
			value[i] = mol[i] / Den;
		}
		return value;
	}

	/**
	 * @param U
	 *            二维数组，指针传递；返回的值
	 * @param V
	 *            一维数组
	 * @param X
	 *            一维数组
	 * @param D
	 *            int
	 * @param uni_array
	 *            二维数组
	 * @param Q
	 *            int
	 * @param J
	 *            int
	 * @return Uniform_local_search List<int[]>
	 */
	public List<int[]> Uniform_local_search(double[][] U, double[] V, double[] X, int D, int[][] uni_array, int Q,
											int J) {
		int[] nouse = new int[D];
		List<int[]> mouse = new ArrayList<int[]>();
		double[][] beta = new double[D][Q];
		zeros(beta);
		zeros(U);

		for (int i = 0; i < D; i++) {
			if (V[i] == X[i]) {				
				ones(beta, i, V[i], false);
			} else {
				double inter = (X[i] - V[i]) / (Q - 1);
				beta[i][0] = V[i];
				beta[i][Q - 1] = X[i];
				for (int k = 1; k < Q - 1; k++) {
					beta[i][k] = beta[i][k - 1] + inter;
				}
			}
		}
		if (D < J) {
			for (int i = 0; i < D; i++) {
				double[] temp = new double[Q];
				temp = beta[i];
				setDoubArray(U, temp, uni_array, i, i);
				mouse.add(i, new int[] { i });
			}
		} else {
			// Divide the variables into J subvectors
			nouse = rnd.randomArray(1, D - 2, D - 2);
			nouse = subArray(nouse, 0, J);
			Arrays.sort(nouse);
			mouse.add(0, getArray(1, nouse[0]));
			for (int i = 1; i < J - 1; i++) {
				mouse.add(i, getArray(nouse[i - 1] + 1, nouse[i]));
			}
			mouse.add(J - 1, getArray(nouse[J - 2] + 1, D));
			// Implement ULS
			for (int i = 0; i < J; i++) {
				double[][] temp = getSubArray(beta, mouse.get(i));
				int k = 0;
				for (int j = 0; j < mouse.get(i).length; j++) {
					double[] temp_1 = temp[j];					
					setDoubArray(U, temp_1, uni_array, i, k);
					k = k + 1;
				}
			}
		}

		return mouse;
	}

	public void setDoubArray(double[][] U, double[] temp, int[][] uni_array, int i, int k) {
		int[] temp1 = new int[uni_array.length];		
		for (int j = 0; j < uni_array.length; j++) {
			temp1[j] = uni_array[j][i];
		}
		double[] temp2 = new double[uni_array.length];
		for (int j = 0; j < temp.length; j++) {
			temp2[j] = temp[temp1[j] - 1];
		}
		for (int j = 0; j < temp2.length; j++) {
			U[j][k] = temp2[j];
		}
	}

	/**
	 * @param beta
	 * @param mouse
	 * @return 
	 */
	public double[][] getSubArray(double[][] beta, int[] mouse) {
		double[][] value = new double[mouse.length][];
		for (int i = 0; i < value.length; i++) {
			value[i] = beta[mouse[i] - 1];
		}
		return value;
	}

	/**
	 * @param i
	 *            开始的值
	 * @param j
	 *            结束的值
	 * @return 返回一个数组 其值为i到j
	 */
	public int[] getArray(int i, int j) {
		int[] value = new int[j - i + 1];
		int line = 0;
		for (int k = i; k < j + 1; k++) {
			value[line] = k;
			line++;
		}
		return value;
	}

	/**
	 * @param value
	 * @param uniform
	 * @param Q 行数 将uniform二维数组中的第Q行去掉
	 */
	public void getSubUnifor(int[][] value, int[][] uniform, int Q) {
		for (int i = 0; i < uniform.length; i++) {
			if (i != Q) {
				for (int j = 0; j < uniform[i].length; j++) {
					value[i][j] = uniform[i][j];
				}
			}
		}
	}

	/**
	 * @param Q
	 *            行的数目
	 * @param J
	 *            列的数目
	 * @return 二维数组
	 */
	public int[][] Uniform_array_generate(int Q, int J) {
		int[][] uni_array = new int[Q][J];
		for (int i = 0; i < Q; i++) {
			uni_array[i][0] = i + 1;
		}
		int i = 2;
		for (int j = 1; j < J; j++) {
			if (gcd(Q, i) != 1) {
				i = i + 1;
			}
			uni_array[0][j] = i;
			i = i + 1;
		}
		for (i = 2; i < Q + 1; i++) {
			for (int j = 2; j < J + 1; j++) {
				if ((i * uni_array[0][j - 1]) % Q == 0) {
					uni_array[i - 1][j - 1] = Q;
				} else {
					uni_array[i - 1][j - 1] = (i * uni_array[0][j - 1]) % Q;
				}
			}
		}

		return uni_array;
	}

	/**
	 * 按轮盘赌策略选择
	 *
	 * @param P
	 *            整数数组
	 * @return 整数1->P.length
	 */
	public int Roulette(double[] P) {
		int Select = 0;
		int m = P.length;
		double r = Math.random();
		double sumP = 0;
		int j = (int) Math.ceil(m * Math.random());
		while (sumP < r) {
			sumP = sumP + P[((j - 1) % m)];
			j = j + 1;
		}
		Select = ((j - 2) % m) + 1;
		return Select;
	}

	/**
	 * @param x
	 * @param y
	 * @return 两个数的公约数
	 */
	public int gcd(int x, int y) {
		if (x == 0)
			return y;
		if (y == 0)
			return x;
		if (x > y)
			return gcd(x % y, y);
		else
			return gcd(x, y % x);
	}

	public void zeros(double[][] array) {
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[i].length; j++) {
				array[i][j] = 0;
			}
		}
	}

	/**
	 * @param array
	 *            数组
	 * @param i
	 *            第几行
	 * @param V
	 *            所赋值
	 * @param t
	 *            布尔类型 false：为行或 ture：为列
	 */
	public void ones(double[][] array, int i, double V, boolean t) {
		for (int j = 0; j < array[i].length; j++) {
			if (t) {
				array[j][i] = V * 1;
			} else {
				array[i][j] = V * 1;
			}
		}
	}

	public void ones(double[][] array) {
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[i].length; j++) {
				array[i][j] = 1;
			}
		}
	}

	public int[] subArray(int[] array, int i, int j) {
		int[] value = new int[j - i];
		for (int k = i; k < j; k++) {
			value[k - i] = array[k];
		}
		return value;
	}

	// 排序
	public void sort(int[] array) {

	}
}
