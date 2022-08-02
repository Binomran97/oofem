package fem;

import inf.v3d.obj.Arrow;
import inf.v3d.obj.Cone;
import inf.v3d.obj.CylinderSet;
import inf.v3d.view.Viewer;

public class Visualizer {
	private double displacementScale;
	private double symbolScale;
	private Structure struct;
	private Viewer viewer;

	public Visualizer(Structure struct, Viewer v) {
		this.struct = struct;
		this.viewer = v;
	}

	public void drawElements() {
		CylinderSet cs = new CylinderSet();

		for (int i = 0; i < this.struct.getNumberOfElements(); i++) {
			double radius = Math.sqrt(this.struct.getElement(i).getArea() / Math.PI);
			double[] p1 = this.struct.getElement(i).getNode1().getPosition().toArray();
			double[] p2 = this.struct.getElement(i).getNode2().getPosition().toArray();
			cs.addCylinder(p1, p2, radius);
		}
		this.viewer.addObject3D(cs);
		// this.viewer.setVisible(true);
	}

	public void drawConstraints() {
		double x;
		double y;
		double z;

		for (int i = 0; i < this.struct.getNumberOfNodes(); i++) {
			Constraint constraint = this.struct.getNode(i).getConstraint();
			x = this.struct.getNode(i).getPosition().getX1();
			y = this.struct.getNode(i).getPosition().getX2();
			z = this.struct.getNode(i).getPosition().getX3();

			for (int j = 0; j < 3; j++) {
				if (constraint.isFree(j) != true) {
					Cone co = new Cone();
					co.setCenter(x, y, z);
					co.setHeight(0.5);
					if (j == 0) {
						co.setDirection(1.0, 0.0, 0.0);
						co.translate(-0.5, 0.0, 0.0);
					}
					if (j == 1) {
						co.setDirection(0.0, 1.0, 0.0);
						co.translate(0.0, -0.5, 0.0);
					}
					if (j == 2) {
						co.setDirection(0.0, 0.0, 1.0);
						co.translate(0.0, 0.0, -0.5);
					}
					this.viewer.addObject3D(co);
				}
			}

		}

		this.viewer.setVisible(true);

	}

	public void drawForces() {
		double x;
		double y;
		double z;
		for (int i = 0; i < this.struct.getNumberOfNodes(); i++) {

			Force forceComponent = this.struct.getNode(i).getForce();
			x = this.struct.getNode(i).getPosition().getX1();
			y = this.struct.getNode(i).getPosition().getX2();
			z = this.struct.getNode(i).getPosition().getX3();
			for (int j = 0; j < 3; j++) {
				if (forceComponent.getComponent(j) != 0) {
					Arrow ar = new Arrow();
					if (j == 0) {
						ar.setPoint1(x - 3.0, y, z);
						ar.setPoint2(x, y, z);
						ar.setRadius(0.05);
					}
					if (j == 1) {
						ar.setPoint1(x, y - 3.0, z);
						ar.setPoint2(x, y, z);
						ar.setRadius(0.05);
					}
					if (j == 2) {
						ar.setPoint1(x, y, z - 3.0);
						ar.setPoint2(x, y, z);
						ar.setRadius(0.05);
					}
					this.viewer.addObject3D(ar);
				}
			}

			this.viewer.setVisible(true);

		}
	}
}
