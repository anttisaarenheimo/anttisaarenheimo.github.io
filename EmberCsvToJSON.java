import java.io.*;
import java.util.*;
import java.util.zip.*;

class EmberCsvToJSON {
	public static void main(String arg[]) throws Exception {
		try {
			if (arg.length != 1) {
				System.out.println("Syntax: 'java EmberCsvToJSON <path to hourly electricity prices zip downloaded from https://ember-energy.org/data/european-wholesale-electricity-price-data/'>");
			}
			ZipFile zip = new ZipFile(arg[0]);
			FileWriter countries = null, options = null;
			for (Enumeration<? extends ZipEntry> entries = zip.entries(); entries.hasMoreElements();) {
				ZipEntry entry = entries.nextElement(); 
				String name = entry.getName();
				if (name.endsWith(".csv") && !name.startsWith("all")) {
					String country = name.substring(0,name.indexOf('.'));
					if (countries == null) {
						options = new FileWriter("options.txt");
						countries = new FileWriter("countries.json");
						countries.write("[\""+country+"\"");
					}
					else {
						countries.write(",\""+country+"\"");
					}
					options.write("<OPTION>"+country+"</OPTION>");
					BufferedWriter json = new BufferedWriter(new FileWriter(name.replace(".csv", ".json")));
					json.write("[");
					BufferedReader reader = new BufferedReader(new InputStreamReader(zip.getInputStream(entry)));
					boolean headerFound = false, isFirst = true;
					String line;
					int	iLocalTime = -1, iPrice = -1;
					while ((line = reader.readLine()) != null) {
						String [] parts = line.split(",");
						if (headerFound && parts.length > iLocalTime && parts.length > iPrice) {
							if (!isFirst) json.write(",");
							json.write("[\""+parts[iLocalTime]+"\","+parts[iPrice]+"]");
							isFirst = false;
						}
						else {
							for (int i = 0; i< parts.length; i++) {
								if (parts[i].indexOf("Local")>=0) iLocalTime = i; 
								else if (parts[i].indexOf("Price")>=0) iPrice = i; 
								headerFound = iLocalTime >= 0 && iPrice >= 0;
								//if (headerFound) System.out.println("Reading "+name+": iLocalTime="+iLocalTime+", iPrice="+iPrice);
							}
						}
					}
					json.write("]");
					json.close();
					reader.close();
				}
			}
			countries.write("]");
			countries.close();
			options.close();
		}
		catch (Exception e) {
			System.out.println(e.toString());
			e.printStackTrace(System.out);
		}
		
	}
}
	