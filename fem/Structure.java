package fem;

import java.util.ArrayList;
import iceb.jnumerics.IMatrix;
import iceb.jnumerics.MatrixFormat;
import iceb.jnumerics.QuadraticMatrixInfo;
import iceb.jnumerics.SolveFailedException;
import iceb.jnumerics.lse.GeneralMatrixLSESolver;
import iceb.jnumerics.lse.ILSESolver;

public class Structure {
	private ArrayList<Element> elements = new ArrayList<Element>();
	private ArrayList<Node> nodes = new ArrayList<Node>();

	public Structure() {

	}

	public Node getNode(int id) {
		return this.nodes.get(id);
	}

	public int getNumberOfNodes() {
		return this.nodes.size();
	}

	public Node addNode(double x1, double x2, double x3) {

		this.nodes.add(new Node(x1, x2, x3));
		return this.nodes.get(getNumberOfNodes() - 1);
	}

	public Element addElement(double e, double a, int n1, int n2) {

		this.elements.add(new Element(e, a, getNode(n1), getNode(n2)));
		return this.elements.get(getNumberOfElements() - 1);
	}

	public Element getElement(int id) {
		return this.elements.get(id);
	}

	public int getNumberOfElements() {
		return this.elements.size();
	}

	private int enumerateDOFs() {
		int NEQ = 0;
		for (int i = 0; i < getNumberOfNodes(); i++) {
			NEQ = getNode(i).enumerateDOFs(NEQ);
		}

		for (int j = 0; j < getNumberOfElements(); j++) {
			getElement(j).enumerateDOFs();
		}
		return NEQ;
	}
	private void assembleStiffnessMatrix(IMatrix kglobalStructure) {
		int [] elementsDOFsNumbers;
		for (int i = 0; i < this.getNumberOfElements(); i++) {
			elementsDOFsNumbers = this.getElement(i).getDOFNumbers();
			for (int j = 0; j < 6; j++) {
				if (elementsDOFsNumbers[j] > -1) {
					for (int k = 0; k < 6; k++) {
						if (elementsDOFsNumbers[k] > -1) {
							kglobalStructure.add(elementsDOFsNumbers[j],elementsDOFsNumbers[k],this.getElement(i).computeElementStiffnessMatrix().get(j, k));
						}
					}
					
				}
			}
		}
		
	}
	private void assembleLoadVector (double[] rGlobalStructure) {
      for (int i = 0; i < this.getNumberOfNodes(); i++) {
    	  for (int j = 0; j < 3 ; j++ ) {
    		  if (getNode(i).getConstraint().isFree(j)) {
    			  rGlobalStructure [getNode(i).getDOFNumbers()[j]] = getNode(i).getForce().getComponent(j);
    		  }
    	  }
      }
	}
	public void solve() {
		int NEQ = this.enumerateDOFs();
		ILSESolver solver = new GeneralMatrixLSESolver();
		QuadraticMatrixInfo aInfo = solver.getAInfo();
		IMatrix a = solver.getA();
		double [] b = new double [NEQ];
		aInfo.setSize(NEQ);
		solver.initialize();
		this.assembleStiffnessMatrix(a);
		this.assembleLoadVector(b);
		
		try {
			solver.solve(b);
		} catch (SolveFailedException e) {
			System.out.println(" Solve failed :  " + e.getMessage());
		}
		
		this.enumerateDOFs();
	}

	public void printResult() {

		System.out.println(this.getNode(getNumberOfNodes() - 1));
	}

}
