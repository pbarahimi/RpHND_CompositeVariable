package RpHND.CompositeVariable.model;
import gurobi.GRBVar;

public class Route implements Comparable<Route>{
	public final Node i;
	public final Node k;
	public final Node m;
	public final Node j;
	public final double cost;
	public final double value;
	public GRBVar var;
	
	/**
	 * Constructor
	 * 
	 * @param h1
	 * @param h2
	 * @param load
	 * @param alpha
	 */
	public Route ( Node i, Node k, Node m, Node j, double alpha ){
		this.i = i;
		this.k = k;
		this.m = m;
		this.j = j;
		this.cost = getRouteCost ( alpha );
		this.value = this.cost * 1 - ( ( 1 - k.failure ) * ( 1 - m.failure ) );
	}
	
	/**
	 * Copy constructor
	 * @param other
	 */
	public Route ( Route other ){
		this.i = other.i;
		this.j = other.j;
		this.k = other.k;
		this.m = other.m;
		this.cost = other.cost;
		this.value = other.value;		
	}
	
	/**
	 * calculates the distance of a route given a load,
	 * origin, destination and a discount factor. 
	 * 
	 * @param load
	 * @param alpha
	 * @return
	 */
	private double getRouteCost( double alpha ){
		double output = getDistance(i, k, alpha);
		output += getDistance(k, m, alpha);
		output += getDistance (m, j, alpha);	
//		System.out.println(load.origin + "_" + hub1 + "_" + hub2 + "_" + load.destination + ": " + output);
		return output;		
	}
	
	/**
	 * returns the Euclidean distance between the two nodes considering
	 * the discount factor if the two nodes are hubs.
	 * @param n1
	 * @param n2
	 * @return
	 */
	private double getDistance (Node n1, Node n2, double alpha ){
		double coefficient = 1;
		if (n1.isHub && n2.isHub) coefficient = 1 - alpha;
		return coefficient * (Math.sqrt(Math.pow(n1.x - n2.x, 2) + Math.pow(n1.y - n2.y, 2)));
	}
	
	@Override
	public String toString(){
		return "(" + this.i + "," + this.k + "," + this.m + "," + this.j + ") - " + this.cost + "-" + this.value ;
	}
	
	@Override
	public int compareTo (Route other){
		if ( this.value < other.value )
			return -1;
		else if ( this.value > other.value )
			return 1;
		else 
			return 0;
	}
}
