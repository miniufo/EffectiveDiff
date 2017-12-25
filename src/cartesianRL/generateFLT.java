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
		float del =3e3f;
		float tstr=93312000;
		float tend=15552000;
		
		List<FltInitData> patch1=generatePatch(new Region2D(2.5e3f,2.5e3f,3075e3f,2195e3f),del,tstr,tend,1000000);
		
		List<FltInitData> all=new ArrayList<>();
		
		all.addAll(patch1);
		
		FLTUtils.toFLTInitFile(all,path+"fltInit_3km_All.bin",DataPrec.float32);
	}
	
	static List<FltInitData> generatePatch(Region2D r,float del,float tstr,float tend,int prefix){
		List<FltInitData> fids=FLTUtils.deployPatch2D(r,del,tstr,tend,1,prefix);
		
		System.out.println("deployed patch ("+prefix+") has "+fids.size()+" particles");
		
		return fids;
	}
}
