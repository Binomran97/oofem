package fem;

import iceb.jnumerics.Vector3D;
import inf.text.ArrayFormat;
import iceb.jnumerics.MatrixFormat;

public class Node {

	private int[] dofNumbers = new int[3];
	private Constraint c;
	private Force f;
	private Vector3D position;
	private Vector3D displacement;

	public Node(double x1, double x2, double x3) {

		this.position = new Vector3D(x1, x2, x3);
        this.c = new Constraint (true, true, true);
        this.f = new Force (0.0, 0.0, 0.0);
	}

	public void setConstraints(Constraint c) {
		this.c = c;

	}

	public Constraint getConstraint() {
		return this.c;

	}

	public void setForce(Force f) {
		this.f = f;
	}

	public Force getForce() {
		return this.f;
	}

	public int[] getDOFNumbers() {
		return this.dofNumbers;
	}

	public int enumerateDOFs(int counter) {
		for (int i = 0; i < 3; i++) {
			Constraint constraint = getConstraint();
			
				if (constraint.isFree(i) == true) {
					this.dofNumbers[i] = counter;
					counter += 1;
				} else {
					this.dofNumbers[i] = -1;
				}
			}
		return counter;
	  }
		
	

	public Vector3D getPosition() {

		return this.position;

	}

//	public void setDisplacement(double u []) {
//		this
//	}
	public void print() {

		System.out.println(MatrixFormat.format(getPosition()));
	}

}
