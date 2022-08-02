package fem;

import java.util.ArrayList;
import iceb.jnumerics.*;
import iceb.jnumerics.lse.*;
import inf.text.ArrayFormat;

public class Structure {

	ArrayList<Node> nodes = new ArrayList<Node>();
	ArrayList<TrussElement> trussElements = new ArrayList<TrussElement>();
	ArrayList<BernoulliElement> bernoulliElements = new ArrayList<BernoulliElement>();
	ArrayList<TimoshenkoElement> timoshenkoElements = new ArrayList<TimoshenkoElement>();

	public Node addNode(double x1, double x2, double x3) {
		Node newNode = new Node(x1, x2, x3);
		nodes.add(newNode);
		return newNode;
	}

	public TrussElement addTrussElement(double e, double a, int n1, int n2) {
		TrussElement newElement = new TrussElement(e, a, nodes.get(n1), nodes.get(n2));
		trussElements.add(newElement);
		return newElement;
	}

	public BernoulliElement addBernoulliElement(double e, double v, double a, int n1, int n2, double iy, double iz) {
		BernoulliElement newElement = new BernoulliElement(e, v, a, nodes.get(n1), nodes.get(n2), iy, iz);
		bernoulliElements.add(newElement);
		return newElement;
	}

	public TimoshenkoElement addTimoshenkoElement(double e, double v, double ky, double kz, double a, int n1, int n2,
			double iy, double iz) {
		TimoshenkoElement newElement = new TimoshenkoElement(e, v, ky, kz, a, nodes.get(n1), nodes.get(n2), iy, iz);
		timoshenkoElements.add(newElement);
		return newElement;
	}

	private void assembleStiffnessMatrix(IMatrix KStructure) {
		int[] elementDoFNumbers;
		for (int i = 0; i < this.getNumberofTrussElements(); i++) {
			elementDoFNumbers = this.getTrussElement(i).getDoFNumbers();
			for (int j = 0; j < 6; j++) {
				if (elementDoFNumbers[j] > -1) {
					for (int k = 0; k < 6; k++) {
						if (elementDoFNumbers[k] > -1) {
							KStructure.add(elementDoFNumbers[j], elementDoFNumbers[k],
									this.getTrussElement(i).computeStiffnessMatrix().get(j, k));
						}
					}
				}
			}
		}
		for (int i = 0; i < this.getNumberofBernoulliElements(); i++) {
			elementDoFNumbers = this.getBernoulliElement(i).getDoFNumbers();
			for (int j = 0; j < 12; j++) {
				if (elementDoFNumbers[j] > -1) {
					for (int k = 0; k < 12; k++) {
						if (elementDoFNumbers[k] > -1) {
							KStructure.add(elementDoFNumbers[j], elementDoFNumbers[k],
									this.getBernoulliElement(i).computeStiffnessMatrix().get(j, k));
						}
					}
				}
			}
		}
		for (int i = 0; i < this.getNumberofTimoshenkoElements(); i++) {
			elementDoFNumbers = this.getTimoshenkoElement(i).getDoFNumbers();
			for (int j = 0; j < 12; j++) {
				if (elementDoFNumbers[j] > -1) {
					for (int k = 0; k < 12; k++) {
						if (elementDoFNumbers[k] > -1) {
							KStructure.add(elementDoFNumbers[j], elementDoFNumbers[k],
									this.getTimoshenkoElement(i).computeStiffnessMatrix().get(j, k));
						}
					}
				}
			}
		}
	}

	private int enumerateDoFs() {
		int neq = 0;
		for (int i = 0; i < nodes.size(); i++) {
			neq = nodes.get(i).enumerateDoFs(neq);
		}
		for (int i = 0; i < bernoulliElements.size(); i++) {
			bernoulliElements.get(i).enumerateDoFs();
		}
		for (int i = 0; i < timoshenkoElements.size(); i++) {
			timoshenkoElements.get(i).enumerateDoFs();
		}
		for (int i = 0; i < trussElements.size(); i++) {
			trussElements.get(i).enumerateDoFs();
		}

		return neq;
	}

	public double getScale() {
		double maxX1 = 0;
		double maxX2 = 0;
		double maxX3 = 0;
		double minX1 = 0;
		double minX2 = 0;
		double minX3 = 0;

		for (int i = 0; i < nodes.size(); i++) {
			if (nodes.get(i).getNode()[0] > maxX1) {
				maxX1 = nodes.get(i).getNode()[0];
			}
			if (nodes.get(i).getNode()[1] > maxX2) {
				maxX2 = nodes.get(i).getNode()[1];
			}
			if (nodes.get(i).getNode()[2] > maxX3) {
				maxX3 = nodes.get(i).getNode()[2];
			}
			if (nodes.get(i).getNode()[0] < minX1) {
				minX1 = nodes.get(i).getNode()[0];
			}
			if (nodes.get(i).getNode()[1] < minX2) {
				minX2 = nodes.get(i).getNode()[1];
			}
			if (nodes.get(i).getNode()[2] < minX3) {
				minX3 = nodes.get(i).getNode()[2];
			}
		}
		double deltaX1 = Math.abs(maxX1 - minX1);
		double deltaX2 = Math.abs(maxX2 - minX2);
		double deltaX3 = Math.abs(maxX3 - minX3);
		double deltaMax = Math.max(deltaX1, Math.max(deltaX2, deltaX3));
		return deltaMax;
	}

	private void assembleLoadVector(double[] rGlobal) {
		for (int i = 0; i < this.getNumberofNodes(); i++) {
			for (int j = 0; j < 6; j++) {
				if (nodes.get(i).getConstraint().isFree(j)) {
					rGlobal[nodes.get(i).getDoFNumbers()[j]] = this.nodes.get(i).getForce().getComponent(j);
				}
			}
		}

	}

	public void solve() {
		int neq = this.enumerateDoFs();
		ILSESolver solver = new GeneralMatrixLSESolver();
		QuadraticMatrixInfo aInfo = solver.getAInfo();
		IMatrix KGlobal = solver.getA();
		double[] rGlobal = new double[neq];
		double[] FGlobal = new double[neq];
		double[] uGlobal = new double[neq];
		aInfo.setSize(neq);
		solver.initialize();
		this.assembleStiffnessMatrix(KGlobal);
		this.assembleLoadVector(rGlobal);
		FGlobal = rGlobal.clone();
		// print
		System.out.println(" Structure Stiffness Matrix");
		System.out.println(MatrixFormat.format(KGlobal));
		System.out.println(" Global Load Matrix");
		System.out.println(ArrayFormat.format(FGlobal));
		try {
			solver.solve(rGlobal);
		} catch (SolveFailedException e) {
			System.out.println(" Solve failed : " + e.getMessage());
		}
		uGlobal = rGlobal.clone();
		Array2DMatrix uGlobalMatrix = new Array2DMatrix(neq, 1, uGlobal);
		Array2DMatrix FGlobalMatrix = new Array2DMatrix(neq, 1, FGlobal);
		IMatrix RGlobalMatrix = (IMatrix) KGlobal.multiply(uGlobalMatrix).subtract(FGlobalMatrix);
		// print result
		System.out.println("  ");
		System.out.println(" Displacement Matrix");
		System.out.println(ArrayFormat.format(uGlobal));
		this.selectDisplacements(rGlobal);
		this.selectNodalForces(RGlobalMatrix);
	}

	private void selectNodalForces(IMatrix RGlobalMatrix) {
		int neq = this.enumerateDoFs();
		double[] RGlobal = new double[neq];
		RGlobal = RGlobalMatrix.getColumn(0).toArray();
		for (int i = 0; i < nodes.size(); i++) {
			double[] nodeForces = { 0, 0, 0, 0, 0, 0 };
			nodes.get(i).getDoFNumbers();
			for (int j = 0; j < 6; j++) {
				if (nodes.get(i).getDoFNumbers()[j] > -1) {
					nodeForces[j] = RGlobal[nodes.get(i).getDoFNumbers()[j]];
				}
			}
			nodes.get(i).setNodalForces(nodeForces);
		}
	}

	private void selectDisplacements(double[] uGlobal) {
		for (int i = 0; i < nodes.size(); i++) {
			double[] nodeDisplacement = { 0, 0, 0, 0, 0, 0 };
			nodes.get(i).getDoFNumbers();
			for (int j = 0; j < 6; j++) {
				if (nodes.get(i).getDoFNumbers()[j] > -1) {
					nodeDisplacement[j] = uGlobal[nodes.get(i).getDoFNumbers()[j]];
				}
			}
			nodes.get(i).setDisplacement(nodeDisplacement);
		}
	}

	public int getNumberofNodes() {
		return nodes.size();
	}

	public Node getNode(int id) {
		return nodes.get(id);
	}

	public int getNumberofTrussElements() {
		return trussElements.size();
	}

	public int getNumberofBernoulliElements() {
		return bernoulliElements.size();
	}

	public int getNumberofTimoshenkoElements() {
		return timoshenkoElements.size();
	}

	public int getNumberofElements() {
		return bernoulliElements.size() + trussElements.size() + timoshenkoElements.size();
	}

	public BernoulliElement getBernoulliElement(int id) {
		return bernoulliElements.get(id);
	}

	public TimoshenkoElement getTimoshenkoElement(int id) {
		return timoshenkoElements.get(id);
	}

	public TrussElement getTrussElement(int id) {
		return trussElements.get(id);
	}

	public void printResults() {
		System.out.println("  ");
		System.out.println("  Listing Structure:");
		System.out.println("  ");
		System.out.println("    Nodes: ");
		System.out.println("       Idx               x1             x2             x3");
		for (int i = 0; i < this.nodes.size(); i++) {
			System.out.print("\t" + i + "\t");
			System.out.println(ArrayFormat.format(nodes.get(i).getPosition().toArray()));

		}
		System.out.println("  ");
		System.out.println("    Constraints: ");
		System.out.println("       Node      u1      u2      u3     phi1    phi2    phi3");
		for (int i = 0; i < this.nodes.size(); i++) {
			if (nodes.get(i).getConstraint().isFree(0) == false
					|| nodes.get(i).getConstraint().isFree(1) == false | nodes.get(i).getConstraint().isFree(2) == false
					|| nodes.get(i).getConstraint().isFree(3) == false
					|| nodes.get(i).getConstraint().isFree(4) == false
					|| nodes.get(i).getConstraint().isFree(5) == false) {
				System.out.print("\t" + i + "\t");
				nodes.get(i).getConstraint().print();
			}
		}
		System.out.println("  ");
		System.out.println("    Forces: ");
		System.out.println(
				"       Node               r1             r2             r3             m1             m2             m3");
		for (int i = 0; i < this.nodes.size(); i++) {
			if (nodes.get(i).getForce().getComponent(0) != 0 || nodes.get(i).getForce().getComponent(1) != 0
					|| nodes.get(i).getForce().getComponent(2) != 0) {
				System.out.print("\t" + i + "\t");
				nodes.get(i).getForce().print();
			}
		}
		System.out.println("  ");
		System.out.println("    Truss Elements: ");
		System.out.println("       Idx               E             A             Length");
		for (int i = 0; i < this.trussElements.size(); i++) {
			System.out.print("\t" + i + "\t");
			trussElements.get(i).print();
		}
		System.out.println("  ");
		System.out.println("    Bernoulli Elements: ");
		System.out.println(
				"       Idx               E             V               A              Iy             Iz           Length");
		for (int i = 0; i < this.bernoulliElements.size(); i++) {
			System.out.print("\t" + i + "\t");
			bernoulliElements.get(i).print();
		}
		System.out.println("  ");
		System.out.println("    Timoshenko Elements: ");
		System.out.println(
				"       Idx               E             V              ky             kz              A              Iy             Iz           Length");
		for (int i = 0; i < this.timoshenkoElements.size(); i++) {
			System.out.print("\t" + i + "\t");
			timoshenkoElements.get(i).print();
		}
		System.out.println("  ");
		System.out.println("  Listing Analysis Results:");
		System.out.println("  ");
		System.out.println("    Displacements: ");
		System.out.println(
				"       Node               u1             u2             u3            ro1            ro2            ro3");
		for (int i = 0; i < this.nodes.size(); i++) {
			System.out.print("\t" + i + "\t");
			System.out.println(ArrayFormat.format(nodes.get(i).getDisplacement()));

		}

		System.out.println("  ");
		System.out.println("    Nodal Forces: ");
		System.out.println(
				"       Node               u1             u2             u3             m1             m2             m3");
		for (int i = 0; i < this.nodes.size(); i++) {
			System.out.print("\t" + i + "\t");
			System.out.println(ArrayFormat.format(nodes.get(i).getNodalForces()));
		}
		System.out.println("  ");
		System.out.println("    Truss Elements Forces: ");
		System.out.println("       Element          Force");
		for (int i = 0; i < this.getNumberofTrussElements(); i++) {
			System.out.print("\t  " + i + "\t");
			System.out.println(trussElements.get(i).computeTrussForce());

		}
	}
}
