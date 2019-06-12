/**
 * COMP3308 Assignment 2, Sem 1 2019
 * kNN and NB algorithms: neither accounts for missing training or testing data
 * @author Ann Cherian
 * SID: 440241928
 */

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class MyClassifier
{
	
	public static int n; // n = number of attributes
	public static ArrayList<double[]> training; // contains all training examples for kNN algorithm
	
	public static ArrayList<ArrayList<Double>> trainingYes; // contains all training examples for NB algorithm where class = yes
	public static ArrayList<ArrayList<Double>> trainingNo; // contains all training examples for NB algorithm where class = no
	public static ArrayList<Double> meanYes; // contains the mean for each attribute with class Yes
	public static ArrayList<Double> meanNo; // contains the mean for each attribute with class No
	public static double Pyes; // P(yes)
	public static double Pno; // P(no)
	public static ArrayList<Double> SDyes; // contains the SD for each attribute with class Yes
	public static ArrayList<Double> SDno; // contains the SD for each attribute with class No
	
	public static void main(String args[]) throws FileNotFoundException
	{		
		// initialise stuff
		n = 0;
		MyClassifier.training = new ArrayList<double[]>();
		MyClassifier.trainingYes = new ArrayList<ArrayList<Double>>();
		MyClassifier.trainingNo = new ArrayList<ArrayList<Double>>();	

		Scanner fil1 = new Scanner(new File(args[0]));		
		
		// find the number of attributes, which is equal to the number of commas in the first line
		String st = fil1.nextLine(); // assume it has at least 1 line
		
		int commas = 0;
		for (int i = 0; i < st.length(); i++)
		{
			if (st.charAt(i) == ',')
			{
				commas++;
			}
		}
		
		n = commas;
				
		Scanner fil2 = new Scanner(new File(args[1]));	
		
		if (args[2].contains("NN"))
		{
			int x = Character.getNumericValue(args[2].charAt(0));
			kNNtraining(st, x, fil1, fil2); // in all cases, the training sample should have at least k examples
		}
		
		else
		{
			//TODO: Run Naive Bayes
			NBtraining(st, fil1, fil2);
		}
	}
	
	/**
	 * Returns the Euclidean distance between two inputs
	 */
	public static double distance(double[] train, double[] test)
	{
		double total = 0;
		
		for (int i = 0; i < n; i++)
		{
			total += Math.pow((double) train[i] - (double) test[i], 2);
		}
		
		return Math.sqrt(total);
	}
	
	/**
	 * Classifies new example based on the nearest k examples in the training data
	 * @param k
	 */
	public static void kNN(int k, Scanner file)
	{		
		while (file.hasNextLine())
		{
			ArrayList<double[]> best = new ArrayList<double[]>(); // contains best k neighbours
			
			// read the test example
			double[] test = new double[n+1];
			String[] str = file.nextLine().split(",");
			
			for (int i = 0; i < n; i++)
			{
				test[i] = Double.parseDouble(str[i]);		
			}
			
			// calculate distance between test and training examples
			for (int i = 0; i < MyClassifier.training.size(); i++)
			{
				double dist = distance(MyClassifier.training.get(i), test);

				if (best.size() < k)
				{
					double[] tmp = {dist, MyClassifier.training.get(i)[n]};
					best.add(tmp);
					// keep best k neighbours sorted, so that the earlier elements are closer neighbours to the test example
					Collections.sort(best, new Comparator<double[]>()
					{
						public int compare(double[] a, double[] b)
						{
							return Double.compare(a[0], b[0]);
						}
					});
				}
				else
				{
					if (dist < best.get(k-1)[0])
					{
						double[] tmp = {dist, MyClassifier.training.get(i)[n]};
						best.set(k-1, tmp);
						// sort best k neighbours. I should create a Comparator function
						Collections.sort(best, new Comparator<double[]>()
						{
							public int compare(double[] a, double[] b)
							{
								return Double.compare(a[0], b[0]);
							}
						});
					}
				}
			}
			// best should contain the best k neighbours
			double vote = 0;
			
			// count number of YES
			for (int i = 0; i < best.size(); i++)
			{
				vote += best.get(i)[1];
			}

			if (vote >= k/2.0)
			{
				System.out.println("yes");
			}
			else
			{
				System.out.println("no");
			}		
		}
	}
	
	/**
	 * Organises the training data for the kNN algorithm
	 * fil1 is for training examples
	 * fil2 is for testing examples (classification)
	 */
	public static void kNNtraining(String st, int k, Scanner fil1, Scanner fil2)
	{
		
		// add the first training example to ArrayList training
		
		double[] first = new double[n+1];
				
		String[] string = st.split(",");

		for (int i = 0; i < n+1; i++)
		{
			if (i == n)
			{
				if (string[i].equals("yes"))
				{
					first[n] = 1.0;
					continue;
				}
				
				else if (string[i].equals("no"))
				{
					first[n] = 0.0;
					continue;
				}				
			}
			else
			{
				first[i] = Double.parseDouble(string[i]);
			}
		}
		MyClassifier.training.add(first);

		while (fil1.hasNextLine())
		{
			double[] temp = new double[n+1];
			
			String[] str = fil1.nextLine().split(",");
			
			for (int i = 0; i < n+1; i++)
			{
				if (i == n)
				{
					if (str[i].equals("yes"))
					{
						temp[n] = 1.0;
						continue;
					}
					
					else if (str[i].equals("no"))
					{
						temp[n] = 0.0;
						continue;
					}
				}
				else
				{
					temp[i] = Double.parseDouble(str[i]);
				}
			}
			MyClassifier.training.add(temp);
		}
		kNN(k, fil2);
	}
	
	public static void NBtraining(String st, Scanner fil1, Scanner fil2)
	{
		// add the first training example to ArrayList training
		
		ArrayList<Double> first = new ArrayList<Double>(n);
				
		String[] string = st.split(",");
		
		if (string[n].equals("yes"))
		{
			for (int i = 0; i < n; i++)
			{
				first.add(Double.parseDouble(string[i]));
			}
			MyClassifier.trainingYes.add(first);
		}
		else
		{
			for (int i = 0; i < n; i++)
			{
				first.add(Double.parseDouble(string[i]));
			}
			MyClassifier.trainingNo.add(first);
		}
		
		while (fil1.hasNextLine())
		{
			ArrayList<Double> temp = new ArrayList<Double>(n);
			
			string = fil1.nextLine().split(",");
			
			if (string[n].equals("yes"))
			{
				for (int i = 0; i < n; i++)
				{
					temp.add(Double.parseDouble(string[i]));
				}
				MyClassifier.trainingYes.add(temp);
			}
			else
			{
				for (int i = 0; i < n; i++)
				{
					temp.add(Double.parseDouble(string[i]));
				}
				MyClassifier.trainingNo.add(temp);
			}
		}
		calculateStats();
		NB(fil2);
	}
	
	/**
	 * Calculates the mean, standard deviation of each attribute, P(yes) and P(no) for NB algorithm
	 */
	public static void calculateStats()
	{
		// calculate mean, SD and P(yes) for all attributes per class
		MyClassifier.meanYes = calcMean(MyClassifier.trainingYes);
		MyClassifier.meanNo = calcMean(MyClassifier.trainingNo);
		
		// calculate P(yes), P(no)
		MyClassifier.Pyes = ((double) MyClassifier.trainingYes.size())/(double) (MyClassifier.trainingYes.size() + MyClassifier.trainingNo.size());
		MyClassifier.Pno = 1 - Pyes;
		
		// calculate SDyes, SDno
		MyClassifier.SDyes = calcSD(MyClassifier.trainingYes, MyClassifier.meanYes);
		MyClassifier.SDno = calcSD(MyClassifier.trainingNo, MyClassifier.meanNo);
	}
	
	public static void NB(Scanner fil2)
	{
		while (fil2.hasNextLine())
		{
			// read the test example
			String[] str = fil2.nextLine().split(",");
			double[] test = new double[n+1];

			for (int i = 0; i < n; i++)
			{
				test[i] = Double.parseDouble(str[i]);		
			}
			
			// TODO:
			double yes = calcCondProb(test, MyClassifier.meanYes, MyClassifier.SDyes, MyClassifier.Pyes);
			double no = calcCondProb(test, MyClassifier.meanNo, MyClassifier.SDno, MyClassifier.Pno);
			
			if (yes >= no)
			{
				System.out.println("yes");
			}
			else
			{
				System.out.println("no");
			}
		}
	}
	
	/**
	 * Calculate the mean for each attribute. The size of ArrayList<Double> = n = the number of attributes
	 * @param arr
	 * @return ArrayList<Double> where each entry is the mean of an attribute
	 */
	public static ArrayList<Double> calcMean(ArrayList<ArrayList<Double>> arr)
	{
		ArrayList<Double> mean = new ArrayList<Double>(n);
		
		// find the sum of each attribute per class
		
		for (int i = 0; i < n; i++)
		{
			mean.add(0.0);
		}
		
		// sum up attribute values
		for (int i = 0; i < arr.size(); i++)
		{
			for (int j = 0; j < n; j++)
			{
				double tmp = mean.get(j);
				mean.set(j, tmp + arr.get(i).get(j));
			}
		}
		
		for (int i = 0; i < n; i++) // loop through each attribute
		{
			mean.set(i, mean.get(i)/((double) arr.size()));			
		}
		
		return mean;
	}
	
	/**
	 * Calculate the SD for each attribute. The size of ArrayList<Double> = n = the number of attributes
	 * @param arr
	 * @param mean
	 * @return ArrayList<Double> where each entry is the SD of an attribute
	 */
	public static ArrayList<Double> calcSD(ArrayList<ArrayList<Double>> arr, ArrayList<Double> mean)
	{
		ArrayList<Double> tmp = new ArrayList<Double>(n); // SD values
		
		if (arr.size() == 1) // only 1 example will mess up SD formula
		{
			for (int i = 0; i < n; i++)
			{
				tmp.add(0.0);
			}
			return tmp;
		}

		for (int i = 0; i < n; i++) // loop through each attribute
		{
			//initialise
			for (int j = 0; j < n; j++)
			{
				tmp.add(0.0);
			}
			
			for (int j = 0; j < arr.size(); j++) // loop through all instances
			{
				double temp = tmp.get(i);
				tmp.set(i, temp + Math.pow(arr.get(j).get(i) - mean.get(i), 2));				
			}
			// we have calculated the sum of (x - mean)^2
		}
		
		for (int j = 0; j < n; j++) // loop through each attribute
		{
			tmp.set(j, Math.sqrt((tmp.get(j)/(double)(arr.size() - 1))));		
		}
		return tmp;
	}
	
	/**
	 * If SD = 0, set conditional probability to 1
	 * @param test = attribute values for the test example
	 * @param mean = meanYES or meanNo
	 * @param SD = SDyes or SDno
	 * @param P = Pyes or Pno
	 * @return P(E1|YES) * P(E2|YES)... * P(E_n|YES) [or ...|NO]
	 */
	public static double calcCondProb(double[] test, ArrayList<Double> mean, ArrayList<Double> SD, double P)
	{
		ArrayList<Double> condProb = new ArrayList<Double>(n);
		// initialise
		for (int i = 0; i < n; i++)
		{
			condProb.add(0.0);
		}
		
		for (int i = 0; i < n; i++)
		{
			if (SD.get(i) == 0.0)
			{
				condProb.set(i, 1.0);
				continue;
			}
			condProb.set(i, Math.exp(-1*(Math.pow(test[i] - mean.get(i), 2))/(2*Math.pow(SD.get(i), 2))));
		}

		for (int i = 0; i < n; i++)
		{
			if (SD.get(i) == 0.0)
			{
				condProb.set(i, 1.0);
				continue;
			}
			condProb.set(i, condProb.get(i)/(SD.get(i)*Math.sqrt(2*Math.PI)));
		}
				
		double numerator = P;
		
		for (int i = 0; i < n; i++)
		{
			numerator *= condProb.get(i);
		}
		return numerator;
	}
}