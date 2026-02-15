# usage: GetEntries() on Reactoribd, alphaN, Geoibd_Th and Geoibd_U 
# ibd: Entries on each run * livetime/(3600*30000), and livetime is calculated from livetime calculator
# (alpha,n): average rate * total livetime
import ROOT
import glob

path = "/data/snoplus3/weiii/antinu/mycuts/Ntuple_afterPosterior/"
ScaledReactoribdpath="/data/snoplus3/weiii/antinu/mycuts/Ntuple_sim/Reactoribd/"
livetimepath = "/data/snoplus3/weiii/antinu/mycuts/Livetime_Calculator/defaultlivetime/"
#eventpath  ={"Reactoribd":f"{ScaledReactoribdpath}scaled_*/","ScaleOccReactoribd":f"{path}scaled_oscillated_ScintFit_2p2Reactoribd/","Geoibd_U":f"{path}Geoibd_U/","Geoibd_Th":f"{path}Geoibd_Th/","Alphan_Lab_13c":f"{path}Alphan_Lab_13c/"}
eventpath={"Reactoribd":f"{ScaledReactoribdpath}scaled_"}
ExpectedNum={"Reactoribd":0.,"ScaleOccReactoribd":0.,"Geoibd_U":0.,"Geoibd_Th":0.,"Alphan_Lab_13c":0.}
day_to_sec = 86400.
def getlivetime(file_path):
	
# Open the file and read the data
	with open(file_path, 'r') as file:
		data = {}
		for line in file:
			key, value = line.split(':')
			data[key.strip()] = float(value.strip())

	# Extract and print the value behind "livetime"
	return float(data["raw_livetime"]), float(data["livetime"])
					


def GetDataFile(listfile):
	RunID = []
	with open(listfile, 'r') as file:
	# Read all lines in the file
		lines = file.readlines()

		# Initialize an empty list to store the first items
		first_items = []

		# Iterate over each line
		for line in lines:
			# Split the line into items based on whitespace
			items = line.split()
			# Append the first item to the list if the line is not empty
			if items:
				RunID.append(items[0])

		# Print the list of first items
		#print(RunID)
	
	return  RunID

if __name__ == "__main__":
	scaling_ibd = 30000.
	Runid = GetDataFile(f"/home/huangp/AntiNu/antinu_runlist_UPDATED.txt")
	for ievtype in eventpath:
		total_livetime = 0.
		total_rawlivetime = 0.
		Tc = ROOT.TChain("DelayT")
		for irun in Runid:
			try:
				raw_livetime,livetime = getlivetime(f"{livetimepath}{irun}livetime_total_defaultlivetime.txt")
				total_livetime += float(livetime)
				total_rawlivetime += float(raw_livetime)
				if ievtype =="Reactoribd":
					Tc.Add(f"{eventpath[ievtype]}{irun}*.root") 
				else:
					Tc.Add(f"{eventpath[ievtype]}*{irun}*.root") 

				
				
				livetime_sec = livetime*day_to_sec
				if(ievtype == "Reactoribd"):
					ExpectedNum[ievtype] += float(Tc.GetEntries())*float(livetime_sec)/(3600.*scaling_ibd)
				else:
					ExpectedNum[ievtype] += float(Tc.GetEntries("post_prob>-3.4"))*float(livetime_sec)/(3600.*scaling_ibd)
				
				Tc.Reset()
			except:
				if(ievtype != "Alphan_Lab_13c"):
					print(f"No matching sim ntuple in run {irun}")
					Tc.Reset()
		
			
		print(f"Expected Num on {ievtype}: {ExpectedNum[ievtype]}")
		print(f"Total livetime on {ievtype} [days]: {total_livetime}")
		print(f"Total raw livetime on {ievtype} [days]: {total_rawlivetime}")
		
