//
package cartesianRL;

import miniufo.util.Region2D;

//
public final class Grids{
	//
	public static final int y=400;
	public static final int x=560;
	
	public static final float res=5500;
	
	public static final float depth=3500;
	public static final float tauMax=0.9f;
	
	public static final String path="D:/Data/MITgcm/barotropicDG/BetaCartRL/";
	
	public static final Region2D region=new Region2D(0,0,res*(x-1),res*(y-1));
}
