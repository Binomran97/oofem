package test;
import fem . Element ;
import fem . Node ;
import iceb . jnumerics .*;
public class StiffnessMatrixTest {
	

	public static void main ( String [] args ) {
	Node n1 = new Node (0, 0, 0);
	Node n2 = new Node (1, 1, 1);
	Element e = new Element (50, Math . sqrt (124) , n1 , n2 );
	IMatrix ke = e. computeElementStiffnessMatrix();
	System .out . println (" Element stiffness matrix ");
	System .out . println ( MatrixFormat . format (ke ));
	}
	}
