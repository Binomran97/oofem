package fem;

	import inf . text . ArrayFormat ;

public class Constraint {

	private boolean[] free = new boolean[3];

	public Constraint(boolean u1, boolean u2, boolean u3) {
		
		this.free[0] = u1;
		this.free[1] = u2;
		this.free[2] = u3;
		
	}
	public boolean isFree(int c) {
	     boolean result;
		 if (this.free[c] == true) {
			 result = true;
		 }
			 else {
				 result = false;
			 }
		 
			 return result;
		 }
		
	public 	String[] isFixedOrfree () {
		String [] state = new String [3];
		for (int i = 0; i < 3; i++) {
			if (this.free[i] == true) {
				state[i] = "Free      ";
			}
			else {
				state [i] = "Fixed     ";
			}
			
		}
		return state;
	}
	
	public void print() {
		System.out.println(ArrayFormat.format(isFixedOrfree()));

	}

}