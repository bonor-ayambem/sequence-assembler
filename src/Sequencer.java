import java.io.*;
import java.util.*;

public class Sequencer {
    public static void main(String[] args){
        List<String> contigs = new ArrayList<String>();

        Scanner sc = new Scanner(System.in);
        List<String> readsList = new ArrayList<String>();

        System.out.println("Enter File Name: ");
        String filename = sc.nextLine();
        File fastaFile = new File(filename);

        // This method call gets the reads from the fasta file
        // and stores them in the list: readsList
        readFile(fastaFile, readsList);
        System.out.println("This file contains " + readsList.size() + " reads");

        // The map variable is the big map that stores all the
        // reads broken up into their prefixes, suffixes, and
        // serial numbers. It will be crucial in assembling the contigs
        Map<Integer, Integer> outdegrees = new HashMap<Integer, Integer>();
        Map<Integer, Integer> indegrees = new HashMap<Integer, Integer>();
        Map<Integer, Map<Integer, List<Integer>>> map = new HashMap<Integer, Map<Integer, List<Integer>>>();
        System.out.println("Creating map from prefixes and suffixes...");
        for(int i = 0; i < readsList.size(); i++){
            int newPrefix = calculateIndex(readsList.get(i).substring(0, 15));
            int newSuffix = calculateIndex(readsList.get(i).substring(15, 30));

            if(map.containsKey(newPrefix)){
                if(map.get(newPrefix).containsKey(newSuffix)){
                    map.get(newPrefix).get(newSuffix).add(i);
                }
                else {
                    map.get(newPrefix).put(newSuffix, new ArrayList<Integer>());
                    map.get(newPrefix).get(newSuffix).add(i);
                }
            }
            else{
                map.put(newPrefix, new HashMap<Integer, List<Integer>>());
                map.get(newPrefix).put(newSuffix, new ArrayList<Integer>());
                map.get(newPrefix).get(newSuffix).add(i);
            }

            if(outdegrees.containsKey(newPrefix)) outdegrees.put(newPrefix, outdegrees.get(newPrefix)+1);
            else outdegrees.put(newPrefix, 1);

            if(indegrees.containsKey(newSuffix)) indegrees.put(newSuffix, indegrees.get(newSuffix)+1);
            else indegrees.put(newSuffix, 1);
        }

        // Finding Paths
        System.out.println("Finding Paths...");
        Set<Integer> usedEdges = new HashSet<Integer>();
        List<Integer> contigNodes = new ArrayList<Integer>();
        boolean pathComplete = true; int startNode = -1;
        while(usedEdges.size() <= readsList.size()){
            if(pathComplete){
                for(Map.Entry<Integer, Map<Integer, List<Integer>>> entry : map.entrySet()){
                    if(!indegrees.containsKey(entry.getKey()) || outdegrees.get(entry.getKey()) > indegrees.get(entry.getKey())){
                        Map<Integer, List<Integer>> subList = map.get(entry.getKey());
                        for(Map.Entry<Integer, List<Integer>> suff : subList.entrySet()){
                            for(int currEdge : suff.getValue()){
                                if(!usedEdges.contains(currEdge)){
                                    startNode = entry.getKey();
                                    pathComplete = false;
                                }
                            }
                        }
                    }
                }
            }

            if(startNode == -1 || !map.containsKey(startNode) || usedEdges.size() == readsList.size()){
                pathComplete = true;
                if(startNode != -1) contigNodes.add(startNode);
                createContig(contigNodes, contigs);
                contigNodes.clear();
                if(usedEdges.size() == readsList.size()) break;
                continue;
            }

            contigNodes.add(startNode);
            int suffixIndex = -1; int unusedEdge = -1;

            boolean found = false;
            Map<Integer, List<Integer>> suffixMap = map.get(startNode);
            for(Map.Entry<Integer, List<Integer>> suffixEntry : suffixMap.entrySet()){
                suffixIndex = suffixEntry.getKey();
                for(int edge : suffixEntry.getValue()){
                    if(!usedEdges.contains(edge)){
                        unusedEdge = edge;
                        found = true;
                        break;
                    }
                }
                if(found) break;
            }

            if(unusedEdge == -1) {
                startNode = -1;
                pathComplete = true;
                createContig(contigNodes, contigs);
                contigNodes.clear();
            }
            else {
                startNode = suffixIndex;
                usedEdges.add(unusedEdge);
            }
        }

        System.out.println("Compressing sequences...");
        while(true){
            List<String> dummy  = new ArrayList<String>(contigs);
            // this method handles instances when one contig is contained in another contig
            // for example, if our result contains one contig: AAAATTTT and another: AAATTT
            // the result should be simply AAAATTTT
            // In other words, if one entire contig is contained inside another one, I will
            // consider it a contig
            accountForDuplicates(contigs);

            // this method handles contigs that overlap
            accountForOffset(contigs);

            Collections.sort(dummy);
            Collections.sort(contigs);

            if(dummy.equals(contigs)) break;
        }

        System.out.println("Writing Results to File...");
        try {
            FileWriter results = new FileWriter("results" + filename);

            for(int i = 0; i < contigs.size(); i++){
                results.write("> contig" + (i+1) + "\n");
                results.write(contigs.get(i) + "\n");
            }
            results.close();
        }
        catch (Exception e) {
            e.getStackTrace();
        }

        System.out.println("PROCESS COMPLETE: Results are located in file " + "results" + filename);
    }

    public static void accountForDuplicates(List<String> contigs){
        Set<String> toRemove = new HashSet<String>();
        Collections.sort(contigs);

        boolean valid = true;
        while(valid) {
            List<String> dummy = new ArrayList<String>(contigs);
            valid = false;
            int i = 0, j = i+1;
            while(i < contigs.size()-1){
                if(areDupes(contigs.get(i), contigs.get(j))){
                valid = true;
                    if(contigs.get(i).length() < contigs.get(j).length())
                        toRemove.add(contigs.get(i));
                    else toRemove.add(contigs.get(j));
                }

                if(j == contigs.size()-1) {
                    i++; j = i+1;
                }
                else j++;
            }

            contigs.removeAll(toRemove);

            Collections.sort(contigs);
            Collections.sort(dummy);
            if(dummy.equals(contigs)) break;
        }
    }
    public static boolean areDupes(String s, String t){
        int i = 0, j = 0;
        List<Integer> nextStart = new ArrayList<Integer>();

        while(i < s.length() && j < t.length()){
            if(!nextStart.contains(i) && s.charAt(i) == t.charAt(0)) nextStart.add(i);
            if(s.charAt(i) == t.charAt(j)){
                if(j == t.length()-1) return true;
                i++; j++;
            }
            else{
                j=0;
                if(nextStart.size() > 1) {
                    i = nextStart.get(1);
                    nextStart.remove(0);
                }
                else i++;
            }
        }

        nextStart.clear();
        // compare s to t
        i = 0; j = 0;
        while(i < t.length() && j < s.length()){
            if(!nextStart.contains(i) && t.charAt(i) == s.charAt(0)) nextStart.add(i);
            if(t.charAt(i) == s.charAt(j)){
                if(j == s.length()-1) return true;
                i++; j++;
            }
            else{
                j=0;
                if(nextStart.size() > 1) {
                    i = nextStart.get(1);
                    nextStart.remove(0);
                }
                else i++;
            }
        }

        return false;
    }

    public static void accountForOffset(List<String> contigs){
        Set<String> toRemove = new HashSet<String>();
        Set<String> toAdd = new HashSet<String>();
        Collections.sort(contigs);

        boolean valid = true;

        while(valid){
            List<String> dummy = new ArrayList<String>(contigs);
            valid = false;
            int i = 0, j = i+1;
            while(i < contigs.size()-1){
                String combo = getOverlap(contigs.get(i), contigs.get(j));
                if(!combo.equals("")){
                valid = true;
                    toRemove.add(contigs.get(i));
                    toRemove.add(contigs.get(j));
                    toAdd.add(combo);
                    i++; continue;
                }

                if(j == contigs.size()-1) {
                    i++; j=i+1;
                }
                else j++;
            }

            contigs.removeAll(toRemove);
            toRemove.clear();

            contigs.addAll(toAdd);
            toAdd.clear();

            Collections.sort(dummy);
            Collections.sort(contigs);

            if(dummy.equals(contigs)) break;
        }
    }
    public static String getOverlap(String s, String t){
        int start = 15;
        int end = Math.min(s.length(), t.length());
        int maxOverlap = 0;

        while(start < end){
            if(s.substring(0, start).equals(t.substring(t.length()-start, t.length())) ||
            s.substring(s.length()-start, s.length()).equals(t.substring(0, start)))
                maxOverlap = start;

            start++;
        }

        if(maxOverlap == 0) return "";

        if(s.substring(0, maxOverlap).equals(t.substring(t.length()-maxOverlap, t.length()))){
            return t + s.substring(maxOverlap, s.length());
        }
        if(s.substring(s.length()-maxOverlap, s.length()).equals(t.substring(0, maxOverlap))){
            return s + t.substring(maxOverlap, t.length());
        }

        return "";
    }

    public static void createContig(List<Integer> nodes, List<String> contigList){
        StringBuilder contig = new StringBuilder();

        for(int i = 0; i < nodes.size(); i++){
            String s = indexToString(nodes.get(i));
            contig.append(s);
        }

        contigList.add(contig.toString());
    }

    public static String indexToString(int index){
        StringBuilder ans = new StringBuilder();

        int power = 14;
        while(power >= 0){
            int val = (int) index / (int)Math.pow(4, power);

            switch(val){
                case 0:
                    ans.insert(0, "A");
                    break;
                case 1:
                    ans.insert(0, "C");
                    break;
                case 2:
                    ans.insert(0, "G");
                    break;
                case 3:
                    ans.insert(0, "T");
                    break;
            }

            index -= val * Math.pow(4, power);
            power--;
        }

        return ans.toString();
    }

    public static int calculateIndex(String read){
        int index = 0;

        for(int i = 0; i < read.length(); i++){
            double curr = 0;

            switch(read.charAt(i)){
                case 'A':
                    curr = 0 * Math.pow(4, i);
                    break;
                case 'C':
                    curr = 1 * Math.pow(4, i);
                    break;
                case 'G':
                    curr = 2 * Math.pow(4, i);
                    break;
                case 'T':
                    curr = 3 * Math.pow(4, i);
                    break;
            }

            index += curr;
        }

        return index;
    }

    public static void readFile(File f, List<String> list){
        try{
            Scanner reader = new Scanner(f);
            while (reader.hasNextLine()) {
                String line = reader.nextLine().trim();
                if (line.charAt(0) != '>') {
                    list.add(line);
                }
            }
        }
        catch(Exception e) {
            System.out.println(e);
        }
    }
}
