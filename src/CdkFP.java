/**
 * Created with IntelliJ IDEA.
 * User: JShen
 * Date: 2/27/13
 * Time: 9:30 PM
 * To change this template use File | Settings | File Templates.
 */
import org.openscience.cdk.Molecule;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.io.*;
import org.openscience.cdk.smiles.*;
import org.openscience.cdk.fingerprint.Fingerprinter;

import java.util.ArrayList;
import java.util.BitSet;
import java.io.*;

public class CdkFP {
    public static ArrayList<String []> readList(String fileName, String separator) throws Exception {
        BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
        String line = null;
        ArrayList<String[]> outList = new ArrayList<String[]>();
        int i=0;
        while((line=br.readLine())!=null){
            String[] data=line.split(separator);
            outList.add(data);
            i++;
        }
        return outList;
    }

    public static void main(String[] args) throws Exception{
        String inputFile="C:\\JShen\\Research\\EDKB\\ER\\DEAC\\130227toGerry\\COMPOUND.txt";
        String outputFile="C:\\JShen\\Research\\EDKB\\ER\\DEAC\\130227toGerry\\COMPOUND_8212_CDKFP.txt";
        ArrayList<String[]> dataMatrix=readList(inputFile,"\t");
        String[] titles=dataMatrix.get(0);
        int smiLoc=0;
        for(int i=0;i<titles.length;i++){
            if (titles[i].equals("SMILES")){
                smiLoc=i;
                break;
            }
        }
        int noInp=dataMatrix.size()-1;

        PrintWriter pw=new PrintWriter(new OutputStreamWriter(new FileOutputStream(outputFile)),true);
        String smi=new String();
        for(int i=1;i<=noInp;i++){
            String[] dataString=dataMatrix.get(i);
            smi=dataMatrix.get(i)[smiLoc];

            for (int j=0;j<titles.length;j++){
                try{
                    pw.printf("%s\t", dataMatrix.get(i)[j]);
                }catch (Exception e){
                    pw.printf("%s\t", "");
                }
            }

            SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
            IAtomContainer mm = sp.parseSmiles(smi);
            Molecule m = new Molecule(mm);
            Fingerprinter fingerprinter = new Fingerprinter();
            BitSet fingerprint = fingerprinter.getFingerprint(m);
            int l=fingerprint.length();
            StringBuffer sb=new StringBuffer (l);

            for (int j=0;j<l;j++){
                if(fingerprint.get(j)){
                    sb.append('1');
                }else{
                    sb.append('0');
                }

            }

            pw.println(sb);
        }
        pw.flush();
        pw.close();
    }



}
