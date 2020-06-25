import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.TreeMap;

import java.util.BitSet;

public class SubFile {
   public static final int BASES_IN=40;
   public static final double THRESHOLD=0.00;
   public static final boolean OUTPUT_OVERLAP = false;
   
   public static final int KEY_ID = 0;

   private class Pair {
      Pair(String id, double beg, double end) {
         this.beg = Math.min(beg,end);
         this.end = Math.max(beg,end);
         this.id = id;
      }
      Pair(String id, double beg, double end, String other) {
         this(id, beg, end);
         info = other;
      }
      public double beg;
      public double end;
      //public double id;
      public String id;
      public String info;
   }
   
   private TreeMap<String, Integer> toRead = new TreeMap<String, Integer>();
   private TreeMap<String, Integer> read = new TreeMap<String, Integer>();
   private LinkedHashMap<String, String> result = new LinkedHashMap<String, String>();
   private TreeMap<String, ArrayList<Pair>> toReadPairs = new TreeMap<String, ArrayList<Pair>>();   
   private TreeMap<String, BitSet> toReadSets = new TreeMap<String, BitSet>();
   
   public double getBasesInSection(String theID, double min, double max) {
      double start = Double.MAX_VALUE;
      double end = Double.MIN_VALUE;
      
      if (true) {
         if (toReadSets.get(theID) == null) {
            ArrayList<Pair> pairs = toReadPairs.get(theID);
            if (pairs == null) { return 0; } // when the id is not meant to be read, skip it

            BitSet b = new BitSet();         
            for (Pair p : pairs) {
               assert(p.id.equalsIgnoreCase(theID));            
               b.set((int)p.beg, (int)p.end);
            }
            toReadSets.put(theID, b);
         }
         BitSet b = toReadSets.get(theID);
         BitSet curr = new BitSet();
         curr.set((int)min, (int)max);
         return (curr.intersects(b) ? 1.0 : -1.0);
         
/*
            double ovl = Math.min(p.end, max) - Math.max(p.beg, min);
            
            if (ovl > 0) {
if (OUTPUT_OVERLAP) {
   NumberFormat nf = new DecimalFormat("#######");
   System.err.println(p.info + " " + nf.format(min-p.beg+1) + " " + nf.format(max-p.beg+1));
}
               if (start > Math.max(p.beg, min)) { start = Math.max(p.beg, min); } 
               if (end < Math.min(p.end, max)) { end = Math.min(p.end, max); }
            }
         }         
*/  
    } else {/*
         for (Double d : toRead.keySet()) {
            if (d < min) {continue;}   
            
            if (d > max) { break; }
            if (start == 0) { start = d; }
            
            if (d <= max) {
               end = d;
            }
         }*/
      }


      return (end-start+1);
   }
   
   public void subFile(String keysFile, String inputFile, Integer idCol, Integer idCol2, boolean rev) throws Exception {
      BufferedReader bf = new BufferedReader(new InputStreamReader(
            new FileInputStream(keysFile)));
      String line = null;
      String prev = null;
      boolean outputNext = false;

System.err.println("BEGIN PROCESSING KEYS");

HashMap<String, ArrayList<String>> lineNum = new HashMap<String, ArrayList<String>>();
      while ((line = bf.readLine()) != null) {
         String[] splitLine = line.trim().split("\\s+");
         
         try {
            if (idCol2 != -1) {               
               //double theID = Double.parseDouble(splitLine[0]);
               String theID = splitLine[0];
               if (toReadPairs.get(theID) == null) {
                  toReadPairs.put(theID, new ArrayList<Pair>());
               }
               toReadPairs.get(theID).add(new Pair(theID, Double.parseDouble(splitLine[1]), Double.parseDouble(splitLine[2]), splitLine[3]));
            } else {
               if (!toRead.containsKey(splitLine[KEY_ID])){
                  toRead.put(splitLine[KEY_ID], 0);  
               }
               //result.put(splitLine[KEY_ID] + toRead.get(splitLine[KEY_ID]).toString(), "");
               toRead.put(splitLine[KEY_ID], toRead.get(splitLine[KEY_ID])+1);               
               if (lineNum.get(splitLine[KEY_ID]) == null) {
                  lineNum.put(splitLine[KEY_ID], new ArrayList<String>());
               }
               lineNum.get(splitLine[KEY_ID]).add(splitLine[0]);
            }
         } catch (Exception e) {
            System.err.println("Warning: Could not parse line " + line + " " + e.getMessage());
         }
      }
System.err.println("DONE PROCESSING KEYS");      
      bf = new BufferedReader(new InputStreamReader(
            new FileInputStream(inputFile)));
      int count = 0;

      while ((line = bf.readLine()) != null) {
         if (line.startsWith("uid")) { System.out.println(line); continue;}
         if (line.startsWith("ID") && !line.startsWith("ID=")) { continue; }
         
         if (count % 1000000 == 0) { System.err.println("PROCESSED " + count + " RECORDS FOR OUTPUT"); }
         String[] splitLine = line.trim().split("\\s+");
         int output = 0;

         if (splitLine.length <= idCol) {
        	 continue;
         }
         String id = splitLine[idCol];//.substring(0, 4);
         if (idCol2 == -1) {
            if (toRead.containsKey(id) || (id.contains("/") && toRead.containsKey(id.substring(0, id.indexOf("/")+1)))) {
               output = toRead.containsKey(id) ? toRead.get(id) : toRead.get(id.substring(0, id.indexOf("/")+1));               
            }
/*
            String altID = id.substring(0, id.length() - 1);
            if (toRead.containsKey(altID)) {
System.err.println("ID " + id + " CHANGED TO " + altID);               
               output = toRead.get(altID);
            }
*/
         }
         else if (idCol2 != -1) {
            double idd = 0;         
            double id2 = -1;
            try {
               idd = Double.parseDouble(splitLine[idCol]);
                  
                  if (idCol2 != -1) {
                     id2 = Double.parseDouble(splitLine[idCol2]);
                  }
            }
            catch (Exception e) {
               System.err.println("WARNING: I could not parse line " + line);
               continue;
            }         
            double min = Math.min(idd,id2);
            double max = Math.max(idd,id2);

            if (true) {
               //Double theID = Double.parseDouble(splitLine[0]);
               String theID = splitLine[splitLine.length-2];
               double bases = getBasesInSection(theID, min, max);
               if ((bases/(max-min+1)) >= THRESHOLD) {
                  output++;                  
               }               
if (output > 0) {
System.err.println("Checking id " + theID + " with " + min + " " + max + " with ovl " + bases + " and output is " + output);
}
            } 
            else {
               double startBases = getBasesInSection("-1", min, min+BASES_IN+1);
               double endBases = getBasesInSection("-1", max-BASES_IN-1, max);
   
               //System.err.println("For range (" + min + ", " + (min+BASES_IN) + ") bases is " + startBases);
               //System.err.println("For range (" + (max-BASES_IN) + ", " + (max) + ") bases is " + endBases);
               if (((startBases/BASES_IN) >= THRESHOLD) && ((endBases/BASES_IN) >= THRESHOLD)) {
                  output++;
               }
            }
         }
         
         if (rev) {
            if (output == 0) { output++; }
            else { output = 0; }
         }

         if (outputNext == true) {
            outputNext = false;
            //System.out.println("THE NEXT " + line);
            //System.out.println("----------------------------------------------------");
         }
         for (int i = 0; i < output; i++) {
            NumberFormat form = NumberFormat.getIntegerInstance();
            form.setGroupingUsed(false);

            //System.out.println("----------- ID " + form.format(id) + " ----------------------");
            //System.out.println("THE PREV " + prev);
            //System.out.println(/*lineNum.get(id).get(i) + " " + */line);
            //result.put(id+i, line);
            System.out.println(line);
            outputNext = true;
            
            Integer val = read.get(id);
            if (val == null) {
               read.put(id, 0);
               val = 0;
            }
            if (!rev) {
               read.put(id, val+1);
            }
            count++;
         }         
         prev = line;
         
         count++;
      }
      
      if (!rev && idCol2 == -1) {
         for (String d : read.keySet()) {
            Integer val = read.get(d);
            
            if (val != null) {
               if (val < 1) {
                  NumberFormat form = NumberFormat.getIntegerInstance();
                  form.setGroupingUsed(false);

                  //System.out.println("ID OF " + form.format(d) + " was read " + val + " times");
                  System.out.println("ID OF " + d + " was read " + val + " times");
               }
            }
         }
      }

      for (String key : result.keySet()) {
         System.out.println(result.get(key));
      }
   }
   
   public static void main(String[] args) throws Exception {
      int idCol = 0;
      int idCol2 = -1;
      boolean rev = false;
      
      if (args.length >= 3) {
         idCol = Integer.parseInt(args[2]);
      }
      if (args.length >=4) {
         idCol2 = Integer.parseInt(args[3]);
      }
      if (args.length >= 5) {
         rev = Boolean.parseBoolean(args[4]);
      }

      SubFile f = new SubFile();
      f.subFile(args[0], args[1], idCol, idCol2, rev);
   }
}
