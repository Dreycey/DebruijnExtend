%%sql
-- your implementation of page rank goes here
--TABLE (nodeid INTEGER, node_group INTEGER, node_symbol VARCHAR(200)) AS
DROP FUNCTION IF EXISTS pagerank();
CREATE OR REPLACE FUNCTION pagerank ()
RETURNS INTEGER AS
$$
DECLARE
    frontier_cursor REFCURSOR;
    select REFCURSOR;
    node_cursor REFCURSOR;
    n_cursor REFCURSOR;
    n_cursor_2 REFCURSOR;
    node_id INTEGER;
    neighbor_id INTEGER;
    check_status INTEGER;
    check_status_1 INTEGER;
    check_status_2 INTEGER;
    connected_component_counter INTEGER;
    symbol_name TEXT;
    neighbor_iteration INTEGER;
    n_neighbor_id INTEGER;
    
BEGIN
 
    --# Create temporary table for all of the nodes.
    DROP TABLE IF EXISTS temp_table;
    CREATE TEMPORARY TABLE temp_table
    AS
    SELECT n.id AS nodeid
    FROM nodes AS n
    LIMIT 0;

    --# Create temporary table for selected nodes. 
    DROP TABLE IF EXISTS selected_table;
    CREATE TEMPORARY TABLE selected_table
    AS
    SELECT n.id AS nodeid, n.id AS node_group, n.symbol AS node_symbol
    FROM nodes AS n
    LIMIT 0;
    
    --# Create temporary table for frontier nodes. 
    DROP TABLE IF EXISTS frontier_table;
    CREATE TEMPORARY TABLE frontier_table
    AS
    SELECT n.id AS nodeid, n.id AS node_iter
    FROM nodes AS n
    LIMIT 0;    
    
    
    OPEN node_cursor FOR EXECUTE 'SELECT n.id FROM nodes as n';
    
    connected_component_counter := 0; --#initialize the counter group ID
    --# Loop for implementing the BFS
    LOOP
        FETCH node_cursor INTO node_id;
        EXIT WHEN NOT FOUND;
        
        --# Check if new node is in selected nodes        
        check_status := (SELECT COUNT(s.nodeid) FROM selected_table AS s WHERE s.nodeid =node_id);

        --# Find maximal connected components if node has not been seen
        IF check_status = 0 THEN
            
            --# Get the symbol name for the node being investigated
            symbol_name := (SELECT n.symbol FROM nodes AS n WHERE n.id=node_id);
            --# Initialize the frontier cursor
            OPEN frontier_cursor FOR EXECUTE format('SELECT e.refid FROM edges as e WHERE e.id=%L',node_id);
            --# add node to the selected table
            --EXECUTE format('INSERT INTO selected_table (nodeid, node_group, node_symbol) VALUES (%L, %L, %L)', node_id, connected_component_counter, symbol_name);
            --# increase the group number for these nodes
            connected_component_counter := connected_component_counter + 1;
            --# Initalize neighbor_iteration
            neighbor_iteration := 0;
            --# Initialize frontier table
            EXECUTE format('INSERT INTO frontier_table (nodeid, node_iter) VALUES (%L, %L)', node_id, neighbor_iteration);
                
            --# Loop through all of the neighbors in the frontier untill empty.
            LOOP
                FETCH frontier_cursor INTO neighbor_id;
                EXIT WHEN NOT FOUND;

                --# Add neighbors reffering to node
                OPEN n_cursor_2 FOR EXECUTE format('SELECT e.id FROM edges as e WHERE e.refid=%L',neighbor_id);
                LOOP
                    FETCH n_cursor_2 INTO n_neighbor_id;
                    EXIT WHEN NOT FOUND;
                    check_status_1 := (SELECT COUNT(s.nodeid) FROM selected_table AS s WHERE s.nodeid =n_neighbor_id);
                    check_status_2 := (SELECT COUNT(f.nodeid) FROM frontier_table AS f WHERE f.nodeid =n_neighbor_id);
                    IF check_status_1 = 0 and check_status_2 = 0 THEN
                        --# ADD TO THE FRONTIER
                        EXECUTE format('INSERT INTO frontier_table (nodeid, node_iter) VALUES (%L, %L)', n_neighbor_id, neighbor_iteration);
                    END IF;
                END LOOP;
                CLOSE n_cursor_2;

                --# Delete current node from frontier temp table
                EXECUTE format('DELETE FROM frontier_table AS f WHERE f.nodeid=%L',neighbor_id);
                --# Update frontier with nodes from updated frontier temp table (ORDER BY neighbor_iteration; most recent first)
                CLOSE frontier_cursor;
                OPEN frontier_cursor FOR EXECUTE 'SELECT f.nodeid FROM frontier_table AS f ORDER BY f.node_iter ASC';
                --# Put the neighbor in the selected, remove from the frontier.
                symbol_name := (SELECT n.symbol FROM nodes AS n WHERE n.id=neighbor_id);
                EXECUTE format('INSERT INTO selected_table (nodeid, node_group, node_symbol) VALUES (%L, %L, %L)', neighbor_id, connected_component_counter, symbol_name);
                
                --# Increase neighbor_iteration by 1 
                neighbor_iteration = neighbor_iteration + 1;
            END LOOP;
            CLOSE frontier_cursor; --# Close the frontier cursor
        END IF; 
    END LOOP;
    CLOSE node_cursor; --# Close the node cursor

    --RETURN QUERY SELECT * FROM selected_table;
    RETURN 0;
END;
$$
LANGUAGE plpgsql;