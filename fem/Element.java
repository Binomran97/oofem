package fem;

import iceb.jnumerics.Array2DMatrix;
import iceb.jnumerics.IMatrix;
import iceb.jnumerics.MatrixFormat;
import iceb.jnumerics.Vector3D;

public class Element {

	private double area;
	private double eModulus;
	private int[] dofNumbers = new int[6];
	private Node node1;
	private Node node2;
	private Vector3D element;

	public Element(double e, double a, Node n1, Node n2) {

		this.area = a;
		this.eModulus = e;
		this.node1 = n1;
		this.node2 = n2;

	}

	public double getArea() {
		return this.area;
	}

	public double getEModulus() {
		return this.eModulus;
	}

	public Node getNode1() {
		return this.node1;
	}

	public Node getNode2() {
		return this.node2;
	}

	public double getLength() {

		double length = (getNode2().getPosition().subtract(getNode1().getPosition())).normTwo();
		return length;
	}

	public Vector3D getE1() {
		this.element = new Vector3D(getEModulus(), getArea(), getLength());
		return element;
	}

	public IMatrix computeElementStiffnessMatrix() {
		double[][] klocalEntries = {
				{ 1 * (this.getEModulus() * getArea() / getLength()), -1 * (getEModulus() * getArea() / getLength()) },
				{ -1 * (getEModulus() * getArea() / getLength()), 1 * (getEModulus() * getArea() / getLength()) } };

		Vector3D positionDifference = getNode2().getPosition().subtract(getNode1().getPosition());
		double[][] transformationMatrixEntries = new double[2][6];

		for (int j = 0; j < 3; j++) {
			transformationMatrixEntries[0][j] = (1 / getLength()) * positionDifference.toArray()[j];
			transformationMatrixEntries[0][j + 3] = 0.0;
			transformationMatrixEntries[1][j] = 0.0;
			transformationMatrixEntries[1][j + 3] = (1 / getLength()) * positionDifference.toArray()[j];
		}

		Array2DMatrix transformationL = new Array2DMatrix(transformationMatrixEntries);
		Array2DMatrix klocal = new Array2DMatrix(klocalEntries);
		Array2DMatrix kglobal = (Array2DMatrix) (transformationL.transpose().multiply(klocal)
				.multiply(transformationL));

		return kglobal;

	}

	public int[] getDOFNumbers() {
		return this.dofNumbers;
	}

	public void enumerateDOFs() {
		for (int i = 0; i < 3; i++) {
			this.dofNumbers[i] = getNode1().getDOFNumbers()[i];
			this.dofNumbers[i + 3] = getNode2().getDOFNumbers()[i];
		}
	}

	public void print() {
		System.out.println(MatrixFormat.format(getE1()));

	}
}
