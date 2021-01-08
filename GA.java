// "CAP'N" Jyym Culpepper
// CMSC 509
// 11/9/10

/* Description:
	
	This program starts by picking 10 random lines of data to use as input
		These random lines come from "parkinsons_updrs.data"
	
	Next, it forms N completely random Equations,
		where N is the number of folds, in this case 10
	Each Equation is formed with 16 Terms
		Each Term is essentially equivalent to a Chromosome for Genetic Algorithm
	
	Next, it begins the Genetic Algorithm, looping MAX_GEN times:
		MAX_GEN is currently set to 5000, as after about 5000 generations,
			any Equation is difficult (and rarely is able) to improve significantly.
		
		Reproduction:
		For every generation, the algorithm reproduces its N "members" X times.
			X is currently set to 4, this can be set to whatever value best suits the algorithm.
			This means that 4 Offspring are generated per Fold.
		
		Crossover:
		For every Term in every Equation, there is a variable probability (C_PROB) to Crossover.
			C_PROB is currently set to 20%.  This is then multiplied by a "factor"
				which limits crossovers as the generations increase.
				This is done so that later generations can avoid crossover when it is not helpful.
			If a Term is selected to Crossover,
				then another Term is selected at random to swap with.
		
		Mutation:
		For every Term in every Equation, there is a fixed probability (M_PROB) to Mutate.
			M_PROB is currently set to 40%.
			If a Term is selected to Mutate,
				then a random Gaussian number is generated to add to the Term.
		
		Selection:
		Every Offspring in every Fold is compared to each other,
			and the best Offspring is selected based on the input data available to this Fold.
	
	Once the Genetic Algorithm has run MAX_GEN times, all Folds are compared to the input data.
		The Average Accuracy is the
			average difference in the best Equation generated in that fold and the actual data.
		
		The Average Accuracy for each fold is displayed.
		The best "Member" is selected.
			The final function of the best "Member" is displayed.
			The Average Accuracy of the best "Member" is displayed.
			The Percent Difference from the actual data is displayed for the best "Member."
*/

import java.io.*;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.util.Random;
import java.text.DecimalFormat;

class GA
{
	public static final int LINES = 5875,   // number of possible inputs
	                        FIELDS = 22,    // number of fields per input
	                        VARIABLES = 16, // number of variables per equation
	                        N = 10,         // number of folds, the N in N-fold
	                        X = 4,          // number of times to replicate an equation
									MAX_GEN = 5000; // number of generations to generate before terminating
	
	public static final double MAX_C = 64,  // maximum (absolute) value a coefficient can start as
	                           MAX_E = 16,  // maximum (absolute) value an exponent can start as
	                           M_PROB = .4, // base probability of Mutation per Term
	                           C_PROB = .2; // base probability of Cross-Over per Term
	
	public static final DecimalFormat D = new DecimalFormat("0.0000"); // For formatting output
	
	public static void main(String[] args)
	{
		int i, j, k, gen = 0, line = 0, newMember, otherFamily, otherOffspring, otherTerm;
		double min, acc=0, prevAcc, avg, factor, newC, newE, totalData = 0.0, data[][] = new double[N][FIELDS];
		Random RNG = new Random();
		
		// pick out N different inputs
		for (i=0;i<N;i++)
		{
			data[i][0] = RNG.nextInt(LINES);
			for (j=0;j<i;j++)
				if (data[j][0] == data[i][0])
				{
					data[i][0] = RNG.nextInt(LINES);
					j = 0;
				}
		}
		sort(data);
		
		// read in said N inputs
		try
		{
			StringTokenizer T;
			//System.out.println("Opening file...");
			Scanner file = new Scanner(new FileReader("parkinsons_updrs.data"));
			System.out.println("Gathering Input...");
			file.nextLine();
			for (i=0;i<N;line++)
			{
				if (data[i][0] != line)
					file.nextLine();
				else
				{
					T = new StringTokenizer(file.nextLine(), ",");
					for (j=0;T.hasMoreTokens();j++)
					{
						data[i][j] = Double.parseDouble(T.nextToken());
						//System.out.print(" " + data[i][j]);
					}
					//System.out.println();
					i++;
				}
			}
			//System.out.println("Closing file...");
			file.close();
			//System.out.println();
		}
		catch(Exception E)
		{
			E.printStackTrace();
			System.exit(0);
		}
		
		System.out.println("Generating Starting Equations...");
		
		// forms N random equations
		double[] coefficient = new double[VARIABLES];
		double[] exponent = new double[VARIABLES];
		Equation[] Member = new Equation[N];
		for (i=0;i<N;i++)
		{
			for (j=0;j<VARIABLES;j++)
			{
				coefficient[j] = 2 * MAX_C * RNG.nextDouble() - MAX_C;
				exponent[j] = 2 * MAX_E * RNG.nextDouble() - MAX_E;
			}
			Member[i] = new Equation(VARIABLES, coefficient, exponent);
		}
		
		// Reproduction, Cross-Over, Mutation, Selection
		System.out.println("Running Genetic Algorithm for " + MAX_GEN + " generations...");
		for(gen=1;gen<=MAX_GEN;gen++)
		{
			// the less factor is, the Members will "update" less, but their updates will be larger
			// the more factor is, the Members will "update" more, but their updates will be smaller
			// in order to overcome those weaknesses, adjust factor as you go
			factor = 50.0 / Math.sqrt(gen); // 100.0 / gen; 
			
			// "Reproduction"
				// reproduce each equation X times, these undergo mutation and crossover
				// reproduce once more, which will be the "baseline" function
			Equation[][] Offspring = new Equation[N][X+1];
			for (k=0;k<N;k++)
				for (i=0;i<X+1;i++)
					Offspring[k][i] = new Equation(Member[k]);
			
			for (k=0;k<N;k++)
			{
				for (i=0;i<X;i++)
				{
					for (j=0;j<VARIABLES;j++)
					{
						// "Cross-Over"
							// may take the term (from this equation) and swap it with the term from another equation
						if (factor * RNG.nextDouble() < C_PROB)
						{
							otherFamily = RNG.nextInt(N);
							otherOffspring = RNG.nextInt(X);
							otherTerm = RNG.nextInt(VARIABLES);
							swap(Offspring[k][i].term[j], Offspring[otherFamily][otherOffspring].term[otherTerm]);
						}
						
						// "Mutation"
							// may mutate the term
						if (RNG.nextDouble() < M_PROB)
						{
							newC = factor * RNG.nextGaussian();
							newE = factor * RNG.nextGaussian();
							Offspring[k][i].term[j].mutate(newC, newE);
						}
					}
				}
			}
			
			// "Selection"
				// Choose equation that best fits input data
			for (k=0;k<N;k++)
			{
				newMember = select(Offspring[k], data, k);
				Member[k] = Offspring[k][newMember];
				
				/*if (gen > 19999)
				{
					acc = totalAccuracy(Member[k], 6, data, k);
					//System.out.println("Generation " + gen + ", Set " + k + " --- Equation " + newMember + ": " + acc);
				}*/
			}
			
			/*newMember = select(Member, data, -1);
			prevAcc = acc;
			acc = totalAccuracy(Member[newMember], 6, data, -1);
			if (prevAcc != acc)
				System.out.println("Generation " + gen + ", best thus far:  Member " + (newMember+1) + " @ " + acc + " total.");*/
		}
		
		System.out.println();
		newMember = 0;
		// min = Member[0].accuracy(6, data[0]);
		min = totalAccuracy(Member[0], 6, data, -1) / 10.0;
		
		// Find the best function
		for (k=0;k<N;k++)
		{
			acc = Member[k].accuracy(6, data[k]);
			avg = totalAccuracy(Member[k], 6, data, -1) / 10.0;
			if (avg < min) // (acc < min)
			{
				min = avg; // = acc;
				newMember = k;
			}
			totalData += data[k][5];
			
			//System.out.println(Member[k]);
			//System.out.println("Accuracy for Member " + (k+1) + ", input " + (k+1) + ": " + D.format(acc));
			System.out.println("Average Accuracy for Member " + (k+1) + ": " + D.format(avg));
		}
		
		System.out.println();
		System.out.println("Best member is Member " + (newMember+1) + ": ");
		System.out.println("\tFunction: " + Member[newMember]);
		// System.out.println("Best member " + (newMember+1) + " accuracy: " + D.format(min));
		System.out.println("\tAverage Accuracy: " + D.format(min));
		System.out.println("\tAverage Percent Difference: " + D.format(1000*min/totalData) + "%");
	}
	
	// sorts random numbers so random lines can be read in properly
	public static void sort(double[][] list)
	{
		int i, j, min, val, temp;
		
		for (i=0;i<N-1;i++)
		{
			min = i;
			val = (int) list[i][0];
			
			for (j=i+1;j<N;j++)
				if (list[j][0] < val)
				{
					min = j;
					val = (int) list[j][0];
				}
			
			if (i != min)
			{
				list[min][0] = list[i][0];
				list[i][0] = val;
			}
		}
	}
	
	// selects best function
	public static int select(Equation[] Offspring, double[][] data, int n)
	{
		int i, newMember = 0;
		double minAcc = totalAccuracy(Offspring[0], 6, data, n), newAcc;
		
		for (i=1;i<Offspring.length;i++)
		{
			newAcc = totalAccuracy(Offspring[i], 6, data, n);
			if (newAcc < minAcc)
			{
				minAcc = newAcc;
				newMember = i;
			}
		}
		
		return newMember;
	}
	
	public static double totalAccuracy(Equation e, int num, double[][] data, int n)
	{
		double total = 0.0;
		
		for (int i=0;i<N;i++)
			if (i != n)
				total += e.accuracy(num, data[i]);
		
		return total;
	}
	
	public static void swap(Term a, Term b)
	{
		Term.swap(a, b);
	}
}

class Equation
{
	private final int VARIABLES;
	public Term[] term;
	
	// Replicates the given equation
	public Equation(Equation Parent)
	{
		this.VARIABLES = Parent.VARIABLES;
		term = new Term[VARIABLES];
		
		for (int i=0;i<VARIABLES;i++)
			this.term[i] = new Term(Parent.term[i]);
	}
	
	// Creates an equation, given a list of terms
	public Equation(int v, Term[] eqn)
	{
		VARIABLES = v;
		term = new Term[v];
		
		for (int i=0;i<v;i++)
			term[i] = eqn[i];
	}
	
	// Creates an equation, given a list of coefficients and exponents
	public Equation(int v, double[] coefficient, double[] exponent)
	{
		VARIABLES = v;
		term = new Term[v];
		
		for (int i=0;i<v;i++)
			term[i] = new Term(coefficient[i], exponent[i]);
	}
	
	// Calculates the equation's accuracy to num (5 is motor_UPDRS, 6 is total_UPDRS)
	public double accuracy(int num, double[] data)
	{
		return Math.abs(data[num-1] - calculate(data));
	}
	
	// Calculates the answer to this equation, with the given data
	public double calculate(double[] data)
	{
		double ans = 0.0;
		
		for (int i=0;i<VARIABLES;i++)
			ans += term[i].calculate(data[data.length-VARIABLES+i]);
		
		return ans;
	}
	
	public String toString()
	{
		String ans = term[0].toString().replace('x', (char) 97);
		
		for (int i=1;i<VARIABLES;i++)
			ans += " + " + term[i].toString().replace('x', (char) (97+i));
		
		return ans;
	}
}

class Term
{
	private double coefficient, exponent;
	private final DecimalFormat D = new DecimalFormat("0.0");
	
	// Replicates the given term
	public Term(Term Parent)
	{
		this.coefficient = Parent.coefficient;
		this.exponent = Parent.exponent;
	}
	
	// Creates a new term, using the given coefficient c and exponent e
	public Term(double c, double e)
	{
		coefficient = c;
		exponent = e;
	}
	
	// Calculates the value of this term (in an equation)
	public double calculate(double data)
	{
		return coefficient * Math.pow(data, exponent);
	}
	
	// Mutates this term's coefficient with c and exponent with e
	public void mutate(double c, double e)
	{
		coefficient += c;
		exponent += e;
	}
	
	// Swaps terms a and b, useful for Cross-Over
	public static void swap(Term a, Term b)
	{
		double c = a.coefficient, e = a.exponent;
		a.coefficient = b.coefficient; a.exponent = b. exponent;
		b.coefficient = c; b.exponent = e;
	}
	
	public String toString()
	{
		return D.format(coefficient) + "*x^" + D.format(exponent);
	}
}
