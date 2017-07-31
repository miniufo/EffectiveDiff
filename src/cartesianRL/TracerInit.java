//
package cartesianRL;

import java.nio.ByteOrder;
import miniufo.basic.ArrayUtil;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.CtlDataWriteStream;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.PI;

//
public final class TracerInit{
	//
	private static final int y=Grids.y;
	private static final int x=Grids.x;
	
	private static final float depth=Grids.depth;
	private static final float tauMax=Grids.tauMax;
	
	private static final String path=Grids.path;
	
	//
	public static void main(String[] args){
		//generateBath(path+"BATH/bath.dat");
		//generateTaux(path+"EXF/taux.dat");
		for(int i=1;i<=16;i++) generateDye(path+"DYE/dye"+i+".dat",i);
	}
	
	static void generateBath(String fname){
		Variable bath=new Variable("bath",new Range(1,1,y,x));
		
		float[][] bdata=bath.getData()[0][0];
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++)
		if(i!=x-1&&j!=y-1) bdata[j][i]=-depth; // set walls
		
		CtlDataWriteStream cdws=new CtlDataWriteStream(fname,ByteOrder.BIG_ENDIAN);
		cdws.writeData(bath); cdws.closeFile();
	}
	
	static void generateDye(String fname,int num){
		Variable dye=new Variable("dye"+num,new Range(1,1,y,x));
		
		float[][] tdata=dye.getData()[0][0];
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++)
		tdata[j][i]=(float)(1.0000312*sin(2.0*PI*j/(y-2.0))*sin(PI*i/(x-2.0))+1.0+num);
		
		System.out.println("dye "+num+": ["+ArrayUtil.getMin(tdata)+", "+ArrayUtil.getMax(tdata)+"]");
		
		CtlDataWriteStream cdws=new CtlDataWriteStream(fname,ByteOrder.BIG_ENDIAN);
		cdws.writeData(dye); cdws.closeFile();
	}
	
	static void generateTaux(String fname){
		Variable taux=new Variable("taux",new Range(1,1,y,x));
		
		float[][] tdata=taux.getData()[0][0];
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++) tdata[j][i]=(float)((1.00001-cos(2.0*PI*j/(y-2.0)))*tauMax/2.0);
		
		CtlDataWriteStream cdws=new CtlDataWriteStream(fname,ByteOrder.BIG_ENDIAN);
		cdws.writeData(taux); cdws.closeFile();
	}
}
