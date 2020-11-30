%Practice file for requesting a file from database within specific date
%range. Useful for seeing what a request looks like, and for practicing
%writing requests while seeing what the error message looks like.
conn = database('Nex_Paradigm_Sort','admin','adminVector123','Vendor','MySql','Server','localhost'); %connect to database
date1 = '2015-11-01';
date2 = '2015-12-01';
subject = 'A';
apolloLMTFselectStatement = ['Select * FROM lmtf WHERE subjID=''', subject,''' AND recDate BETWEEN CAST(''', date1, ''' AS DATE) AND CAST(''', date2, ''' AS DATE);'];
apolloRequest = exec(conn, apolloLMTFselectStatement);
if size(apolloRequest.Message,1)==0
    apolloFetch = fetch(apolloRequest);
    columnNames = columnnames(apolloFetch);
    if size(apolloRequest.Message,1)==0
        ApolloNovemberData = apolloFetch.Data;
    else
        disp(apolloRequest.Message);
    end
else
    disp(apolloRequest.Message);
end
disp(['Columns: ', columnNames]);
disp(ApolloNovemberData);
close(apolloRequest); close(apolloFetch); close(conn);