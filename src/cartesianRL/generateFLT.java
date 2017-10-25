//
package cartesianRL;

import java.util.ArrayList;
import java.util.List;
import common.MITgcmUtil.DataPrec;
import fltDispersion.FLTUtils;
import fltDispersion.FltInitData;
import miniufo.util.Region2D;

//
public final class generateFLT{
	//
	private static final String path=Grids.path;
	
	
	//
	public static void main(String[] args){
		float del =2e3f;
		float tstr=31104000;
		float tend=93312000;
		
		List<FltInitData> patch1=generatePatch(new Region2D( 400e3f, 1050e3f, 404e3f, 1150e3f), del, tstr, tend, 10000);
		List<FltInitData> patch2=generatePatch(new Region2D(2500e3f, 1050e3f,2504e3f, 1150e3f), del, tstr, tend, 20000);
		List<FltInitData> patch3=generatePatch(new Region2D(2500e3f, 1900e3f,2504e3f, 2000e3f), del, tstr, tend, 30000);
		List<FltInitData> patch4=generatePatch(new Region2D(2500e3f,  400e3f,2504e3f,  500e3f), del, tstr, tend, 40000);
		List<FltInitData> patch5=generatePatch(new Region2D(1600e3f, 1450e3f,1604e3f, 1550e3f), del, tstr, tend, 50000);
		List<FltInitData> patch6=generatePatch(new Region2D(1500e3f,  100e3f,1504e3f,  200e3f), del, tstr, tend, 60000);
		
		List<FltInitData> all=new ArrayList<>();
		
		all.addAll(patch1);
		all.addAll(patch2);
		all.addAll(patch3);
		all.addAll(patch4);
		all.addAll(patch5);
		all.addAll(patch6);
		
		FLTUtils.toFLTInitFile(all,path+"flt_init.bin",DataPrec.float32);
	}
	
	static List<FltInitData> generatePatch(Region2D r,float del,float tstr,float tend,int prefix){
		List<FltInitData> fids=FLTUtils.deployPatch2D(r,del,tstr,tend,1,prefix);
		
		System.out.println("deployed patch ("+prefix+") has "+fids.size()+" particles");
		
		return fids;
	}
}
