/* 
Dreycey Albin
Declarative SQL
functions
for COMP5333
project: DIPSS 
*/







--This query was for finding proteins involved with apoptosis.
--NOTE: This will not work in the most recent DIPSS database,
--as the go tables have been removed from the UniProt table
--to normalize the database.

CREATE VIEW joinjoin AS select uni.*, h.pdb_id
FROM uniprot as uni INNER JOIN hasproteins as h ON uni.uniprot_id=h.uniprot_id;
--where publicationyear > 2017;

select ss.pdb_id, SUM(ss.proteinlength)
FROM joinjoin AS jj INNER JOIN secondarystructure as ss ON jj.pdb_id=ss.pdb_id
--WHERE jj.proteinlength > 300
WHERE jj.goprocess LIKE '%apoptosis%'
GROUP BY ss.pdb_id;

SELECT count(*)
FROM UniProt;







--Gathering the number of proteins with length 40, and the number of chains

/*
This can be useful to the user as one may be interested in collecting proteins
above (or below) a certain length.
This can then be grouped by the PDB, giving back the number of chains.
Depending on the length, change from 40. 
*/

SELECT pdb_id, ROUND(SUM(proteinlength)/40) AS numberOfChains
FROM SecondaryStructure
WHERE proteinlength = 40
GROUP BY pdb_id;






--Asking the quesiton: What is residuelength?
SELECT CS.pdb_id, residuecount, sum(proteinlength)
FROM crystalStructureSequence AS CS 
INNER JOIN SECONDARYSTRUCTURE AS SS 
ON CS.pdb_id = SS.pdb_id 
GROUP BY CS.pdb_id, residuecount
ORDER BY residuecount DESC;






--Looking at joins with publicaiton year included. 
SELECT COUNT(*) 
FROM CrystalStructure CS
INNER JOIN SECONDARYSTRUCTURE AS SS 
ON CS.pdb_id = SS.pdb_id 
INNER JOIN crystalStructureSequence AS CSS
ON CSS.pdb_id = CS.pdb_id 
WHERE publicationyear < 1980 AND NOT publicationyear=-1






--What go terms are associated with most proteins?
SELECT gofunctionid, COUNT(gofunctionid) AS numberofproteins
FROM GOFUNCTIONS
group by gofunctionid
ORDER BY numberofproteins DESC;





--How man unique GO terms are in the table?
SELECT count(*) 
FROM (SELECT gofunctionid
		FROM GOFUNCTIONS 
		group by gofunctionid) AS foo;




		
--Select for proteins involved with apoptosis
SELECT *
FROM GOPROCESSES AS GP
WHERE GP.CELLPROCESS LIKE '%apoptosis%';






--How man unique GO processes are in the table? --> 3331
SELECT count(*) 
FROM (SELECT cellprocessid
		FROM GOPROCESSES
		group by cellprocessid) AS foo;

		
SELECT kmer, GP.uniprot_id, SS.PDB_ID, CHAINCODE, SST3
FROM GOPROCESSES AS GP
INNER JOIN UNIPROT AS U 
ON U.uniprot_id = GP.uniprot_id
INNER JOIN HasProteins AS HP
ON HP.uniprot_id = U.uniprot_id
INNER JOIN secondaryStructure AS SS
ON SS.pdb_id = HP.pdb_id;

 
SELECT GP.uniprot_id, SS.PDB_ID, CHAINCODE, SST3, SS.PROTEINSEQUENCE, GP.cellProcessid, GP.cellProcess
FROM GOPROCESSES AS GP
INNER JOIN UNIPROT AS U 
ON U.uniprot_id = GP.uniprot_id
INNER JOIN HasProteins AS HP
ON HP.uniprot_id = U.uniprot_id
INNER JOIN secondaryStructure AS SS
ON SS.pdb_id = HP.pdb_id
WHERE GP.uniprot_id LIKE '%A0A024RAV5%' --GP.CELLPROCESS LIKE '%apoptosis%'
GROUP BY GP.uniprot_id, SS.PDB_ID, Chaincode, SST3, SS.PROTEINSEQUENCE, GP.cellProcessid, GP.cellProcess
 



