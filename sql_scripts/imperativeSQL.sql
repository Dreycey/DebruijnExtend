/* 
Dreycey Albin
Imperative SQL
functions
for COMP5333
project: DIPSS 
*/






--############################################
--Imperative SQL for grabbing kmers from a input string.
--############################################

DROP FUNCTION get_kmers(TEXT, integer);

CREATE OR REPLACE FUNCTION get_kmers (proteinSequence TEXT, kmerlength INT) 
RETURNS TABLE (k_mers TEXT) 
AS 
$$
DECLARE
   counter INTEGER := 1 ; 
   sequenceLength INT;
   
BEGIN

--SELECT LEN('AGTGCTAGCTAGC');
SELECT char_length(proteinSequence) INTO sequenceLength;
sequenceLength := sequenceLength - (kmerlength - 1);

LOOP

RETURN QUERY
SELECT SUBSTRING(proteinSequence, counter, kmerlength);
EXIT WHEN counter = sequenceLength;
counter := counter + 1;
EXIT WHEN kmerlength > sequenceLength;
END LOOP; 

END; 
$$ 
LANGUAGE 'plpgsql';

--############ END OF FUNCTION #####################

select * from get_kmers('AGGCAGTCGATGCTGATCGATGACTGATGCATCGATCGATCGATCGTAGATCGTAGTAGCTGACTGACTCG', 900)

--Using the following view for tett runs
CREATE VIEW testuni AS SELECT * 
FROM UNIPROT
LIMIT 10;






--############################################
--Imperative SQL for finding proteins with similiar kmers 
--############################################

DROP FUNCTION findlikeproteins(proteinSequence TEXT, kmerlength INT);

CREATE OR REPLACE FUNCTION findlikeproteins(proteinSequence TEXT, kmerlength INT) 
RETURNS TABLE (k_mers TEXT, matchingprotein VARCHAR(3000)) 
AS 
$$
DECLARE
   sequenceLength INT;
   queryString TEXT;
   counter INTEGER := 1;
   kmerTable CURSOR FOR
   SELECT * FROM 
   get_kmers(proteinSequence, kmerlength);
   kmer TEXT;
BEGIN
SELECT char_length(proteinSequence) INTO sequenceLength;
sequenceLength := sequenceLength - (kmerlength - 1);

OPEN kmerTable;
LOOP
        FETCH kmerTable INTO kmer; EXIT WHEN NOT FOUND;
        RETURN QUERY
        SELECT kmer, uniprot_id FROM TESTUNI AS U WHERE U.proteinsequence LIKE '%'||kmer||'%';
        counter := counter + 1;
END LOOP;
CLOSE kmerTable;
END; 
$$ 
LANGUAGE 'plpgsql';

--############ END OF FUNCTION #####################

--Using the function

select * from findlikeproteins('MLAFIPVLTKKMNPRSTEAAIKY', 21);

--Interesting query to runL
select matchingprotein from findlikeproteins('MLAFIPVLTKKMNPRSTEAAIKY', 21)
GROUP BY matchingprotein;





--############################################
--Imperative SQL for grabbing potential molecular functions from kmers of a input string. 
--############################################

DROP FUNCTION findPotentialMolecularFunctions(proteinSequence TEXT, kmerlength INT);

CREATE OR REPLACE FUNCTION findPotentialMolecularFunctions(proteinSequence TEXT, kmerlength INT) 
RETURNS TABLE (k_mers TEXT, GO_id VARCHAR(3000), MolecularFunction VARCHAR(3000)) 
AS 
$$
DECLARE
   sequenceLength INT;
   queryString TEXT;
   counter INTEGER := 1;
   kmerTable CURSOR FOR
   SELECT * FROM 
   get_kmers(proteinSequence, kmerlength);
   kmer TEXT;
BEGIN
SELECT char_length(proteinSequence) INTO sequenceLength;
sequenceLength := sequenceLength - (kmerlength - 1);

OPEN kmerTable;
LOOP
        FETCH kmerTable INTO kmer; EXIT WHEN NOT FOUND;
        RETURN QUERY
        SELECT kmer, GF.gofunctionid, GF.gofunction FROM UNIPROT AS U 
        INNER JOIN GOFUNCTIONS AS GF ON U.uniprot_id = GF.uniprot_id
        WHERE U.proteinsequence LIKE '%'||kmer||'%';
        counter := counter + 1;
END LOOP;
CLOSE kmerTable;
END; 
$$ 
LANGUAGE 'plpgsql';

--############ END OF FUNCTION #####################


--An interesting query to run with this function:
select go_id, MolecularFunction from findPotentialMolecularFunctions('MLAFIPVAIKY', 5)
GROUP BY go_id, MolecularFunction;

--WOAH..
SELECT go_id, MolecularFunction, COUNT(go_id) AS Hits 
FROM findPotentialMolecularFunctions('MDEADRRLLRRCRLRLVEELQVDQLWDALLSRELFRPHMIEDIQRAGSGSRRDQARQLII
DLETRGSQALPLFISCLEDTGQDMLASFLRTNRQAAKLSKPTLEN', 20)
GROUP BY go_id, MolecularFunction
ORDER BY Hits DESC
LIMIT 3;

--Testing a combination of proteins, showing the power of using Kmers.
SELECT go_id, molecularFunction, COUNT(go_id) AS Hits 
FROM findPotentialMolecularFunctions('MNENLFASFIAPTILGLPAAVCRLRLVEELQVDQLWDALLSRELFRPHMIEDIQRAGSGSRRDQARQLII
DLETRGSQALPLFISCLEDTGQDMLASFLRTNRQALIILFPPLLIPTSK', 5)
GROUP BY go_id, MolecularFunction
ORDER BY Hits DESC
LIMIT 10;







--############################################
--Imperative SQL for grabbing potential cell functions from kmers of a input string. 
--############################################

DROP FUNCTION findPotentialCellularFunctions(proteinSequence TEXT, kmerlength INT);

CREATE OR REPLACE FUNCTION findPotentialCellularFunctions(proteinSequence TEXT, kmerlength INT) 
RETURNS TABLE (k_mers TEXT, GO_id VARCHAR(3000), CellularFunction VARCHAR(3000)) 
AS 
$$
DECLARE
   sequenceLength INT;
   queryString TEXT;
   counter INTEGER := 1;
   kmerTable CURSOR FOR
   SELECT * FROM 
   get_kmers(proteinSequence, kmerlength);
   kmer TEXT;
BEGIN
SELECT char_length(proteinSequence) INTO sequenceLength;
sequenceLength := sequenceLength - (kmerlength - 1);

OPEN kmerTable;
LOOP
        FETCH kmerTable INTO kmer; EXIT WHEN NOT FOUND;
        RETURN QUERY
        SELECT kmer, GP.cellprocessid, GP.cellprocess FROM UNIPROT AS U 
        INNER JOIN GOPROCESSES AS GP ON U.uniprot_id = GP.uniprot_id
        WHERE U.proteinsequence LIKE '%'||kmer||'%';
        counter := counter + 1;
END LOOP;
CLOSE kmerTable;
END; 
$$ 
LANGUAGE 'plpgsql';

--############ END OF FUNCTION #####################


--An interesting query to run with this function:
select go_id, cellularFunction from findPotentialCellularFunctions('MLAFIPVAIKY', 5)
GROUP BY go_id, cellularFunction;

--WOAH..
SELECT go_id, cellularFunction, COUNT(go_id) AS Hits 
FROM findPotentialCellularFunctions('MDEADRRLLRRCRLRLVEELQVDQLWDALLSRELFRPHMIEDIQRAGSGSRRDQARQLII
DLETRGSQALPLFISCLEDTGQDMLASFLRTNRQAAKLSKPTLEN', 20)
GROUP BY go_id, cellularFunction
ORDER BY Hits DESC
LIMIT 3;

--Testing a combination of proteins, showing the power of using Kmers.
SELECT go_id, cellularFunction, COUNT(go_id) AS Hits 
FROM findPotentialCellularFunctions('MNENLFASFIAPTILGLPAAVCRLRLVEELQVDQLWDALLSRELFRPHMIEDIQRAGSGSRRDQARQLII
DLETRGSQALPLFISCLEDTGQDMLASFLRTNRQALIILFPPLLIPTSK', 5)
GROUP BY go_id, cellularFunction
ORDER BY Hits DESC
LIMIT 10;


--############################################
--Wrapper needed to find the index of a substring within a bigger string:
--############################################

CREATE OR REPLACE FUNCTION "patindex"( "pattern" VARCHAR, "expression" VARCHAR ) RETURNS INT AS $BODY$
SELECT
    COALESCE(
        STRPOS(
             $2
            ,(
                SELECT
                    ( REGEXP_MATCHES(
                        $2
                        ,'(' || REPLACE( REPLACE( TRIM( $1, '%' ), '%', '.*?' ), '_', '.' ) || ')'
                        ,'i'
                    ) )[ 1 ]
                LIMIT 1
            )
        )
        ,0
    )
;
$BODY$ LANGUAGE 'sql' IMMUTABLE;

--############ END OF FUNCTION #####################

--FROM: https://stackoverflow.com/questions/43914647/sql-patindex-equivalent-in-postgresql
--Similiar to the PATINDEX here: https://www.w3schools.com/sql/func_sqlserver_patindex.asp

--Imperative SQL for grabbing secondary structure possibilites from kmers of a input string. 
SELECT PATINDEX('%Sool%', 'W3Schools.com');

--SELECT SUBSTRING('SQL Tutorial', PATINDEX('%hool%', 'W3Schools.com'), 8) AS ExtractString;




--############################################
-- Function for finding an input proteins secondary structure.
--############################################

DROP FUNCTION findsimiliarsecondarystructures(text,integer);

CREATE OR REPLACE FUNCTION findSimiliarSecondaryStructures(proteinSequence TEXT, kmerlength INT) 
RETURNS TABLE (k_mers TEXT, indexx INT, Guniprot_id VARCHAR(3000), proteinsequence VARCHAR(3000), PDB_ID VARCHAR(3000), CHAINCODE VARCHAR(3000), SST3 TEXT, SST8 TEXT) 
AS 
$$
DECLARE
   sequenceLength INT;
   queryString TEXT;
   counter INTEGER := 1;
   kmerTable CURSOR FOR
   SELECT * FROM 
   get_kmers_2(proteinSequence, kmerlength);
   kmer TEXT;
BEGIN
SELECT char_length(proteinSequence) INTO sequenceLength;
sequenceLength := sequenceLength - (kmerlength - 1);

OPEN kmerTable;
LOOP
        FETCH kmerTable INTO kmer, indexx; EXIT WHEN NOT FOUND;
        RETURN QUERY
                SELECT kmer, indexx, GP.uniprot_id, GP.proteinsequence, SS.PDB_ID, SS.CHAINCODE, 
                                SUBSTRING(SS.SST3, PATINDEX(kmer, U.proteinsequence), kmerlength), 
                                SUBSTRING(SS.SST8, PATINDEX(kmer, U.proteinsequence), kmerlength)
                FROM GOPROCESSES AS GP
                INNER JOIN UNIPROT AS U 
                ON U.uniprot_id = GP.uniprot_id
                INNER JOIN HasProteins AS HP
                ON HP.uniprot_id = U.uniprot_id
                INNER JOIN secondaryStructure AS SS
                ON SS.pdb_id = HP.pdb_id
                WHERE U.proteinsequence LIKE '%'||kmer||'%'
                ORDER BY indexx;
        counter := counter + 1;
END LOOP;
CLOSE kmerTable;
END; 
$$ 
LANGUAGE 'plpgsql';

--INteresting query to run:
SELECT indexx, sst3, count(sst3)
from findSimiliarSecondaryStructures('GSHMKNSVSVDLPGSMKVLVSKSSNADGKYDLIATVDALELSGTSDKNNGSGVLEGVKADASKVKLTISD', 5)
--WHERE sst3 not like '' AND length(sst3) = 5
GROUP BY indexx, sst3
ORDER BY indexx, count(sst3) DESC


--############################################
--Imperative SQL for grabbing kmers from a input string. V2
--Used in a different way, with different output, for predicting
--the proteins secondary structures. 
--############################################

DROP FUNCTION get_kmers_2(TEXT, integer);

CREATE OR REPLACE FUNCTION get_kmers_2 (proteinSequence TEXT, kmerlength INT) 
RETURNS TABLE (k_mers TEXT, indexx INT) 
AS 
$$
DECLARE
   counter INTEGER := 1 ; 
   sequenceLength INT;
   
BEGIN

--SELECT LEN('AGTGCTAGCTAGC');
SELECT char_length(proteinSequence) INTO sequenceLength;
sequenceLength := sequenceLength - (kmerlength - 1);

LOOP

RETURN QUERY
SELECT SUBSTRING(proteinSequence, counter, kmerlength), counter;
EXIT WHEN counter = sequenceLength;
counter := counter + 1;
EXIT WHEN kmerlength > sequenceLength;
END LOOP; 

END; 
$$ 
LANGUAGE 'plpgsql';

--############ END OF FUNCTION #####################



select * from get_kmers_2('AGGCAGTCGATGCTGATCGATGACTGATGCATCGATCGATCGATCGTAGATCGTAGTAGCTGACTGACTCG', 900)



--############################################
--Function/algorithm for predicting possible secondary structures:
--############################################

DROP FUNCTION predictSecondaryStructure(proteinSequence TEXT, kmerlength INT);

CREATE OR REPLACE FUNCTION predictSecondaryStructure(proteinSequence TEXT, kmerlength INT) 
RETURNS TABLE (proteinSecondarySequence TEXT, weightTotal BIGINT, proteinlength INT) AS

$$
DECLARE
   sequenceLength INT;
   queryString TEXT;
   counter INTEGER := 1;
   createviewquery TEXT;
   CURSOR_1 REFCURSOR;
   CURSOR_2 REFCURSOR;
   CURSOR_3 REFCURSOR;
   counterLength INT;
   createtemporarytablequery TEXT;
   querystring_1 TEXT;
   querystring_2 TEXT;
   queryString_3 TEXT;
   proteinSecondarySequence TEXT;
   dropViewQuery TEXT;
   dropTableQuery TEXT;
   weightTotal BIGINT;
   proteinlength INT;
   indexx INT;
   indexx_2 INT;
   sst3struct TEXT; 
   queryMatches BIGINT;
   extensionTableLength INT;
   insertQuery_1 TEXT;
   insertQuery TEXT;
   indexx_last INT DEFAULT 0;
   endOfGrowingStruct VARCHAR(3000);
   startOfStruct VARCHAR(3000);
   indexStart INT;
   counterLengthQuery TEXT;
   endOfStruct TEXT;
   newproteinSecondarySequence VARCHAR(3000);
   newweightTotal BIGINT;
   
BEGIN
SELECT char_length(proteinSequence) INTO sequenceLength;
sequenceLength := sequenceLength - (kmerlength - 1);

--Create a view specifically for this function to work with as the input table

RAISE NOTICE 'NOW WE IS Starting' ; 
DROP TABLE IF EXISTS secondaryStructView;
CREATE TEMPORARY TABLE secondaryStructView
AS
SELECT SS.indexx, SS.sst3, count(SS.sst3)
from findSimiliarSecondaryStructures( QUOTE_LITERAL(proteinSequence), kmerlength) AS SS
GROUP BY SS.indexx, SS.sst3
ORDER BY SS.indexx, count(SS.sst3) DESC;

EXECUTE 'SELECT indexx FROM secondaryStructView LIMIT 1' INTO indexStart;


--This temporary table will hold the potential sequences as they are being
--extended by the algorithm
DROP TABLE IF EXISTS extensionTable;
CREATE TEMPORARY TABLE extensionTable 
AS
SELECT SS.sst3 AS proteinSecondarySequence,
SS.count AS weightTotal,
CHAR_LENGTH(SS.sst3) AS proteinlength
FROM secondaryStructView AS SS
WHERE SS.indexx = indexStart;
--(proteinSecondarySequence TEXT, weightTotal INT, proteinlength INT);

counterLength := kmerlength; --Set initial length of the predicted
                                                
--First loop goes through the different indexes:
queryString_3 = 'SELECT indexx FROM secondaryStructView 
WHERE indexx NOT IN (SELECT indexx
FROM secondaryStructView 
ORDER BY indexx
LIMIT 1)
GROUP BY indexx 
ORDER BY indexx';

OPEN CURSOR_3 FOR EXECUTE queryString_3;
                                                                                                                  
LOOP
        FETCH CURSOR_3 INTO indexx; EXIT WHEN NOT FOUND;
                                                                                                                  
        --Getting cursor for first nested loop prepared:
        queryString_1 = 'SELECT *
                         FROM secondaryStructView
                         WHERE indexx =' || indexx;
        OPEN CURSOR_1 FOR EXECUTE queryString_1;
        RAISE NOTICE 'NOW WE IS AT SECOND LOOP' ;
                RAISE NOTICE 'The index is % ', indexx;
        --This will loop through the kmers associated with index from first loop:                                                          
        LOOP
                FETCH CURSOR_1 INTO indexx_2, sst3struct, queryMatches; EXIT WHEN NOT FOUND;
                
                
                /*This will loop will ensure that the k-mer matches the end 
                of the sequence being extended. If there's a match, it will
                add the secondary structure notation to the end. If not, it
                will go back through the loop.*/
                queryString_2 = 'SELECT *
                                 FROM extensionTable
                                 WHERE proteinlength = ' || counterLength;             
                OPEN CURSOR_2 FOR EXECUTE queryString_2;
                                                                                                                  
                LOOP
                        EXECUTE 'SELECT COUNT(*) FROM extensionTable' INTO extensionTableLength;
                        RAISE NOTICE 'BEFORE: COUNTERLENGTH IS AT %, %, % ', counterLength, extensionTableLength, kmerlength;
                        FETCH CURSOR_2 INTO proteinSecondarySequence, weightTotal, proteinlength; EXIT WHEN NOT FOUND;
                                                
                                                --Making variables to make more readable
                                                endOfGrowingStruct = SUBSTRING(proteinSecondarySequence, 
                                                                                   (CHAR_LENGTH(proteinSecondarySequence)-(kmerlength-2)),(kmerlength-1));    
                                                startOfStruct = SUBSTRING(sst3struct, 1, (kmerlength-1));
                                endOfStruct = SUBSTRING(sst3struct, (CHAR_LENGTH(sst3struct) - (indexx-indexx_last-1)), (indexx-indexx_last)); 
                                                
                                                RAISE NOTICE 'BEFORE: SEQUENCES BEING COMPARED ARE %, % ', 
                                                        endOfGrowingStruct, startOfStruct;
                                                        
                                                newweightTotal = weightTotal + queryMatches;
                                                newproteinSecondarySequence = proteinSecondarySequence || endOfStruct;
                                                 IF endOfGrowingStruct = startOfStruct THEN
                                insertQuery = format('INSERT INTO extensionTable(proteinSecondarySequence, weightTotal, proteinlength) 
                                               VALUES (%L, %L, %L)',
                                               newproteinSecondarySequence, newweightTotal, CHAR_LENGTH(newproteinSecondarySequence));
                                EXECUTE insertQuery;
                         END IF;
                        RAISE NOTICE 'AFTER: COUNTERLENGTH IS AT %', counterLength;
                                                
                        
                END LOOP;
                CLOSE CURSOR_2;
        END LOOP;
        CLOSE CURSOR_1;
                counterLengthQuery = 'SELECT MAX(ET.proteinlength) FROM extensionTable AS ET';
        EXECUTE counterLengthQuery INTO counterLength;
                indexx_last := indexx;
END LOOP;
CLOSE CURSOR_3;

--CLOSE kmerTable;
RETURN QUERY SELECT * FROM extensionTable;
--Drop the view for the function
RAISE NOTICE 'AT SECOND DROP #1 %', counterLength;
DROP TABLE IF EXISTS secondaryStructView;
--Drop the temporary table for the function
RAISE NOTICE 'AT SECOND DROP #2 %', counterLength;
DROP TABLE IF EXISTS extensionTable; 
END; 
$$ 
LANGUAGE 'plpgsql';

--############ END OF FUNCTION #####################



SELECT * FROM predictSecondaryStructure('MLAFIPVLTKKMNPRSTEAAIKY', 4);







--Build this later:
--select * from predictSecondaryStructureGoBoost('MLAFIPVLTKKMNPRSTEAAIKY', 21);







