% The sheet we are working on has to be moved to the first position in
% order to be able to read the data. This does not apply to writing to
% sheets.

% uMontreal analysis
clientId    = '946432324902-5l38jed11bdrvmq0e8emu4f9cspfs9en.apps.googleusercontent.com';
clientSecret= 'KOb1uzbFNMMkUz5pT6NHzOgV';
RunOnce(clientId, clientSecret);

%%
uMontrealDataDir = 'C:\Users\Mark\Dropbox\Projet EIT';
addpath(genpath(uMontrealDataDir));
URL1 = 'https://docs.google.com/spreadsheets/d/18Pgu7nbQPOQSEUqBum-U-paW0ZLBsDvfr4ankbQAJmA/edit#gid=0';
URL2 = 'https://docs.google.com/spreadsheets/d/18Pgu7nbQPOQSEUqBum-U-paW0ZLBsDvfr4ankbQAJmA/edit#gid=1333875953';
URL3 = 'https://docs.google.com/spreadsheets/d/18Pgu7nbQPOQSEUqBum-U-paW0ZLBsDvfr4ankbQAJmA/edit#gid=449260216';
status = start_from_scratch(URL2);

%%
% modify breath boundaries
[str,stat]  = urlread('https://docs.google.com/spreadsheets/d/18Pgu7nbQPOQSEUqBum-U-paW0ZLBsDvfr4ankbQAJmA/export?format=csv');
data        = web_string_to_table(str);
newTriads   = parse_triads(data);
if isempty(newTriads)
    eitFile = load_eit(fileName);
    status  = populate_sheet(spreadsheetID, sheetID, eitFile);
else
    eitFile.triads = newTriads;
end

show_breath_boundaries(eitFile.fdata, triads);
%%

function [data, status] = get_sheet_data(spreadsheetID, sheetName, sheetpos, d)
    % include java classes
    import java.io.*;
    import java.net.*;
    import java.lang.*;
    com.mathworks.mlwidgets.html.HTMLPrefs.setProxySettings
    aSheets = refreshAccessToken; % refresh and retrieve access token
    % parameters
    url = ['https://sheets.googleapis.com/v4/spreadsheets/' spreadsheetID '/values/' sheetName];
    maxiter = 3;
    status = 1;

    str_d = savejson('',d); % builds values string in JSON
    str_range = BuildRange(sheetpos, d); % given data d and starting pos, builds range 

    if size(d,1)==1, str_d = ['[' str_d ']']; end

    success = 0;
    iter = 1;
    while ~success && iter<maxiter    
        iter=iter+1;
        con = urlreadwrite(mfilename,url);
        con.setInstanceFollowRedirects(false);
        con.setRequestMethod('GET');
        con.setDoOutput(true);
        con.setDoInput(true);
        con.setRequestProperty('Authorization',['Bearer ' aSheets]);
        con.setRequestProperty('Content-Type','application/json; charset=UTF-8');               
        con.setRequestProperty('X-Upload-Content-Length', '0');
        event = ['{',...
                    '''valueInputOption'': ''RAW'',',...
                    '''data'': [{',...
                        '''range'': ''' sheetName '!' str_range ''',',...
                        '''values'': ' str_d ',',...
                    '}],',...
                 '}'];

        con.setRequestProperty('Content-Length', num2str(length(event)));

        ps = PrintStream(con.getOutputStream());
        ps.print(event);
        ps.close();  clear ps;  
        if (con.getResponseCode()~=200)
            con.disconnect();
            continue;
        end
        success=true;
    end

    if ~success
        status = 0;
        display(['Failed pasting data. Last response was: ' num2str(con.getResponseCode) '/' con.getResponseMessage().toCharArray()']);
    end
end


function status = start_from_scratch(URL)
    [spreadsheetID, sheetID] = parse_sheet_url(URL);
    aSheets = refreshAccessToken; % refresh and retrieve access token
    [~, ~, sheetName] = getSheetMetaData(spreadsheetID, sheetID, aSheets); % to ensure we have enough cells
    fileName = horzcat(sheetName, '.eit');
    eitFile = load_eit(fileName);
    status = populate_sheet(spreadsheetID, sheetID, eitFile);
%     update_image();
end


function [spreadsheetID, sheetID] = parse_sheet_url(URL)
    splitURL    = strsplit(URL, '/');
    spreadsheetID = splitURL{end-1};
    sheetID     = strsplit(splitURL{end}, 'gid=');
    sheetID     = sheetID{end};
end


function plot_ventilation(eitFile)
    global fig_path;
%     global axsFS;
    global titleFS;
    fig     = bigFig(1);
    triads  = eitFile.triads;
    nBreaths= size(triads,1);
    if nBreaths > 0
        nRows = 4;
        subplot( nRows, nBreaths, 1:nBreaths*nRows-nBreaths);
        show_breath_boundaries(eitFile.fdata, triads);
        calc_colours('cmap_type', 'blue_red');
        colormap(calc_colours('colourmap'));
        hold on;
        for j=1:nBreaths
            subplot( nRows, nBreaths, nBreaths*nRows-nBreaths+j);
            image(calc_breath_delta_z(imgs, triads(j,:)).*1000 + 128);
            title(num2str(j));
            axis equal; axis off;
        end
        hold off;
    else
        show_breath_boundaries(eitFile.fdata, triads);
    end
    sgtitle(remove_underscores(eitFile.name), 'FontSize',titleFS);
    printPDF( sprintf('%s%s Global Z', fig_path, eitFile.name) );
end


function show_breath_boundaries(data, triads)
    % show seconds
    global textFS;
    temp = triads;
    gsig = sum(data,1);
    textY = min(gsig);
    plot(gsig); hold on;
    for i = 1:size(temp,1)
       plot(temp(i,:), gsig(temp(i,:)), 'or');
       text(temp(i,2), textY, num2str(i), 'FontSize',textFS);
    end
end


function status = populate_sheet(spreadsheetID, sheetID, eitFile)
    triads  = horse_breath_finder(eitFile.fdata);
    status = mat2sheets(spreadsheetID, sheetID, [1,1], {'EELI1','EILI','EELI2'});
    status = mat2sheets(spreadsheetID, sheetID, [2,1], triads);
end


function EIT = load_eit(file, startIdx, endIdx)
    if nargin < 3
        endIdx = -1;
    end
    if nargin ==1
        startIdx = 1;
    end
    
    [vv,aux]= eidors_readdata(file);
    if endIdx == -1
        endIdx = size(vv,2);
    end
    t_rel   = aux.t_rel(:, startIdx:endIdx);
    fs      = 1e6 ./ median( diff(t_rel) ); %framerate is median dif of time points/ 1000000 (convert to s)
    vv      = vv(:, startIdx:endIdx);
    clip    = 50;
    rvv     = real(vv);
    fdata   = lowpass(rvv',1, fs,'StopbandAttenuation',40,'ImpulseResponse','iir')';
    fdata   = fdata(:, clip: (size(fdata,2) - clip) ); % trim filter edge artifacts
    EIT     = struct;
    EIT.fs  = fs;
    EIT.data    = vv;
    EIT.fdata   = fdata;
end % end function


function triads = parse_triads(webData)
    if isempty(webData)
        triads = [];
        return
    end
    triads = nan(size(webData,1)-1, 3);
    header = webData(1,:);
    for i=1:length(header)
        if strcmp(header(i), 'EELI1')
            triads(:,1)  = str2double(webData(2:end,i));
        elseif strcmp(header(i), 'EILI')
            triads(:,2)  = str2double(webData(2:end,i));
        elseif strcmp(header(i), 'EELI2')
            triads(:,3)  = str2double(webData(2:end,i));
       end
    end
end


function data = web_string_to_table(str)
    if isempty(str)
        data = [];
        return
    end
    data = textscan(str,'%s','whitespace','\n');
    data = data{1};
    for ii=1:length(data)
       %for each line, split the string into its comma-delimited units
       %the '%q' format deals with the "quoting" convention appropriately.
       tmp = textscan(data{ii},'%q','delimiter',',');
       data(ii,1:length(tmp{1})) = tmp{1};
    end
end


function aSheets = refreshAccessToken
    aSheets = [];

    if ~exist('google_tokens.mat', 'file')
        disp('Run RunOnce() first and ensure that google_tokens.mat is in MATLAB''s path'); 
        return
    end

    load google_tokens.mat

    newAccessTokenString=urlread('https://accounts.google.com/o/oauth2/token','POST', ...
    {'client_id', client_id, 'client_secret', client_secret, 'refresh_token', rSheets, 'grant_type', 'refresh_token'});

    aSheets=[];

    reply_commas=[1 strfind(newAccessTokenString,',') length(newAccessTokenString)];

    for i=1:length(reply_commas)-1
        if ~isempty(strfind(newAccessTokenString(reply_commas(i):reply_commas(i+1)),'access_token'))
            tmp=newAccessTokenString(reply_commas(i):reply_commas(i+1));
            index_tmp_colon=strfind(tmp,':');
            tmp=tmp(index_tmp_colon+1:end); clear index_tmp_colon;
            index_quotes=find(tmp=='"');
            aSheets=tmp(index_quotes(1)+1:index_quotes(2)-1); clear index_quotes tmp;
        end
    end

    save('google_tokens.mat', 'aSheets', '-append');
end

function [nRows, nCols, sheetName] = getSheetMetaData(spreadsheetID, sheetID, aSheets)
    % include java classes
    import java.io.*;
    import java.net.*;
    import java.lang.*;
    com.mathworks.mlwidgets.html.HTMLPrefs.setProxySettings

    % params
    %url = ['https://spreadsheets.google.com/feeds/worksheets/' sheetID '/private/full/' spreadsheetID];
    url = ['https://sheets.googleapis.com/v4/spreadsheets/' spreadsheetID '?includeGridData=false'];
    MAXITER=10;
    success=false;
    safeguard=0;

    %init
    nRows = [];
    nCols = [];
    sheetName = [];

    while (~success && safeguard<MAXITER)
        safeguard=safeguard+1;
        con = urlreadwrite('mat2sheets',url);
        con.setInstanceFollowRedirects(false);
        con.setRequestMethod('GET');
        con.setDoInput(true);
        con.setRequestProperty('Authorization',['Bearer ' aSheets]);       

        if (con.getResponseCode()~=200)
            con.disconnect();
            continue;
        end
        success=true;
    end

    if success
        json_return = [];
        isr = java.io.InputStreamReader(con.getInputStream);
        br = java.io.BufferedReader(isr);
        l = br.readLine();
        while ~isempty(l)
            json_return = [json_return, l];
            l = br.readLine();
        end
        j = char(json_return)';
        j = j(:)';

        opt.SimplifyCell = 1;
        jdat = loadjson(j, opt);
        sIds = nan(size(jdat.sheets));
        if iscell(jdat.sheets(1)) % test for some weird bug
            for i = 1:numel(jdat.sheets) 
                jd = jdat.sheets(i); 
                jdc = cell2mat(jd); 
                sIds(i) = jdc.properties.sheetId; 
            end 

            wheresheet = find(sIds==str2double(sheetID));

            if ~isempty(wheresheet) && length(wheresheet)==1 
                jd = jdat.sheets(wheresheet); 
                jdc = cell2mat(jd); 
                nRows = jdc.properties.gridProperties.rowCount; 
                nCols = jdc.properties.gridProperties.columnCount; 
                sheetName = jdc.properties.title; 
            else 
                disp('SheetID doesnt match sheets'); 
            end 
        else
            for i = 1:numel(jdat.sheets)
                sIds(i) = jdat.sheets(i).properties.sheetId;
            end

            wheresheet = find(sIds==str2double(sheetID));

            if ~isempty(wheresheet) && length(wheresheet)==1
                nRows = jdat.sheets(wheresheet).properties.gridProperties.rowCount;
                nCols = jdat.sheets(wheresheet).properties.gridProperties.columnCount;
                sheetName = jdat.sheets(wheresheet).properties.title;
            else
                disp('SheetID doesnt match sheets');
            end
        end

        con.disconnect();

    else
        display(['Failed gathering metadata. Last response was: ' num2str(con.getResponseCode) '/' con.getResponseMessage().toCharArray()']);
    end
end


function [urlConnection,errorid,errormsg] = urlreadwrite(fcn,urlChar)
    %URLREADWRITE A helper function for URLREAD and URLWRITE.

    %   Matthew J. Simoneau, June 2005
    %   Copyright 1984-2009 The MathWorks, Inc.
    %   $Revision: 1.1.6.5 $ $Date: 2009/09/28 20:28:02 $

    % Default output arguments.
    urlConnection = [];
    errorid = '';
    errormsg = '';

    % Determine the protocol (before the ":").
    protocol = urlChar(1:min(find(urlChar==':'))-1);

    % Try to use the native handler, not the ice.* classes.
    switch protocol
        case 'http'
            try
                handler = sun.net.www.protocol.http.Handler;
            catch exception %#ok
                handler = [];
            end
        case 'https'
            try
                handler = sun.net.www.protocol.https.Handler;
            catch exception %#ok
                handler = [];
            end
        otherwise
            handler = [];
    end

    % Create the URL object.
    try
        if isempty(handler)
            url = java.net.URL(urlChar);
        else
            url = java.net.URL([],urlChar,handler);
        end
    catch exception %#ok
        errorid = ['MATLAB:' fcn ':InvalidUrl'];
        errormsg = 'Either this URL could not be parsed or the protocol is not supported.';
        return
    end

    % Get the proxy information using MathWorks facilities for unified proxy
    % preference settings.
    mwtcp = com.mathworks.net.transport.MWTransportClientPropertiesFactory.create();
    proxy = mwtcp.getProxy(); 


    % Open a connection to the URL.
    if isempty(proxy)
        urlConnection = url.openConnection;
    else
        urlConnection = url.openConnection(proxy);
    end
end


function data = loadjson(fname,varargin)
    %
    % data=loadjson(fname,opt)
    %    or
    % data=loadjson(fname,'param1',value1,'param2',value2,...)
    %
    % parse a JSON (JavaScript Object Notation) file or string
    %
    % authors:Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
    % created on 2011/09/09, including previous works from 
    %
    %         Nedialko Krouchev: http://www.mathworks.com/matlabcentral/fileexchange/25713
    %            created on 2009/11/02
    %         François Glineur: http://www.mathworks.com/matlabcentral/fileexchange/23393
    %            created on  2009/03/22
    %         Joel Feenstra:
    %         http://www.mathworks.com/matlabcentral/fileexchange/20565
    %            created on 2008/07/03
    %
    % $Id$
    %
    % input:
    %      fname: input file name, if fname contains "{}" or "[]", fname
    %             will be interpreted as a JSON string
    %      opt: a struct to store parsing options, opt can be replaced by 
    %           a list of ('param',value) pairs - the param string is equivallent
    %           to a field in opt. opt can have the following 
    %           fields (first in [.|.] is the default)
    %
    %           opt.SimplifyCell [0|1]: if set to 1, loadjson will call cell2mat
    %                         for each element of the JSON data, and group 
    %                         arrays based on the cell2mat rules.
    %           opt.FastArrayParser [1|0 or integer]: if set to 1, use a
    %                         speed-optimized array parser when loading an 
    %                         array object. The fast array parser may 
    %                         collapse block arrays into a single large
    %                         array similar to rules defined in cell2mat; 0 to 
    %                         use a legacy parser; if set to a larger-than-1
    %                         value, this option will specify the minimum
    %                         dimension to enable the fast array parser. For
    %                         example, if the input is a 3D array, setting
    %                         FastArrayParser to 1 will return a 3D array;
    %                         setting to 2 will return a cell array of 2D
    %                         arrays; setting to 3 will return to a 2D cell
    %                         array of 1D vectors; setting to 4 will return a
    %                         3D cell array.
    %           opt.ShowProgress [0|1]: if set to 1, loadjson displays a progress bar.
    %
    % output:
    %      dat: a cell array, where {...} blocks are converted into cell arrays,
    %           and [...] are converted to arrays
    %
    % examples:
    %      dat=loadjson('{"obj":{"string":"value","array":[1,2,3]}}')
    %      dat=loadjson(['examples' filesep 'example1.json'])
    %      dat=loadjson(['examples' filesep 'example1.json'],'SimplifyCell',1)
    %
    % license:
    %     BSD License, see LICENSE_BSD.txt files for details 
    %
    % -- this function is part of JSONLab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
    %

    global pos inStr len  esc index_esc len_esc isoct arraytoken

    if(regexp(fname,'^\s*(?:\[.+\])|(?:\{.+\})\s*$','once'))
       string=fname;
    elseif(exist(fname,'file'))
       try
           string = fileread(fname);
       catch
           try
               string = urlread(['file://',fname]);
           catch
               string = urlread(['file://',fullfile(pwd,fname)]);
           end
       end
    else
       error('input file does not exist');
    end

    pos = 1; len = length(string); inStr = string;
    isoct=exist('OCTAVE_VERSION','builtin');
    arraytoken=find(inStr=='[' | inStr==']' | inStr=='"');
    jstr=regexprep(inStr,'\\\\','  ');
    escquote=regexp(jstr,'\\"');
    arraytoken=sort([arraytoken escquote]);

    % String delimiters and escape chars identified to improve speed:
    esc = find(inStr=='"' | inStr=='\' ); % comparable to: regexp(inStr, '["\\]');
    index_esc = 1; len_esc = length(esc);

    opt=varargin2struct(varargin{:});

    if(jsonopt('ShowProgress',0,opt)==1)
        opt.progressbar_=waitbar(0,'loading ...');
    end
    jsoncount=1;
    while pos <= len
        switch(next_char)
            case '{'
                data{jsoncount} = parse_object(opt);
            case '['
                data{jsoncount} = parse_array(opt);
            otherwise
                error_pos('Outer level structure must be an object or an array');
        end
        jsoncount=jsoncount+1;
    end % while

    jsoncount=length(data);
    if(jsoncount==1 && iscell(data))
        data=data{1};
    end

    if(isfield(opt,'progressbar_'))
        close(opt.progressbar_);
    end
end


function opt=varargin2struct(varargin)
    %
    % opt=varargin2struct('param1',value1,'param2',value2,...)
    %   or
    % opt=varargin2struct(...,optstruct,...)
    %
    % convert a series of input parameters into a structure
    %
    % authors:Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
    % date: 2012/12/22
    %
    % input:
    %      'param', value: the input parameters should be pairs of a string and a value
    %       optstruct: if a parameter is a struct, the fields will be merged to the output struct
    %
    % output:
    %      opt: a struct where opt.param1=value1, opt.param2=value2 ...
    %
    % license:
    %     BSD License, see LICENSE_BSD.txt files for details 
    %
    % -- this function is part of jsonlab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
    %

    len=length(varargin);
    opt=struct;
    if(len==0) return; end
    i=1;
    while(i<=len)
        if(isstruct(varargin{i}))
            opt=mergestruct(opt,varargin{i});
        elseif(ischar(varargin{i}) && i<len)
            opt=setfield(opt,lower(varargin{i}),varargin{i+1});
            i=i+1;
        else
            error('input must be in the form of ...,''name'',value,... pairs or structs');
        end
        i=i+1;
    end
end


function s=mergestruct(s1,s2)
    %
    % s=mergestruct(s1,s2)
    %
    % merge two struct objects into one
    %
    % authors:Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
    % date: 2012/12/22
    %
    % input:
    %      s1,s2: a struct object, s1 and s2 can not be arrays
    %
    % output:
    %      s: the merged struct object. fields in s1 and s2 will be combined in s.
    %
    % license:
    %     BSD License, see LICENSE_BSD.txt files for details 
    %
    % -- this function is part of jsonlab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
    %

    if(~isstruct(s1) || ~isstruct(s2))
        error('input parameters contain non-struct');
    end
    if(length(s1)>1 || length(s2)>1)
        error('can not merge struct arrays');
    end
    fn=fieldnames(s2);
    s=s1;
    for i=1:length(fn)              
        s=setfield(s,fn{i},getfield(s2,fn{i}));
    end
end


function val=jsonopt(key,default,varargin)
    %
    % val=jsonopt(key,default,optstruct)
    %
    % setting options based on a struct. The struct can be produced
    % by varargin2struct from a list of 'param','value' pairs
    %
    % authors:Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
    %
    % $Id: loadjson.m 371 2012-06-20 12:43:06Z fangq $
    %
    % input:
    %      key: a string with which one look up a value from a struct
    %      default: if the key does not exist, return default
    %      optstruct: a struct where each sub-field is a key 
    %
    % output:
    %      val: if key exists, val=optstruct.key; otherwise val=default
    %
    % license:
    %     BSD License, see LICENSE_BSD.txt files for details
    %
    % -- this function is part of jsonlab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
    % 

    val=default;
    if(nargin<=2) return; end
    opt=varargin{1};
    if(isstruct(opt))
        if(isfield(opt,key))
           val=getfield(opt,key);
        elseif(isfield(opt,lower(key)))
           val=getfield(opt,lower(key));
        end
    end
end


function parse_char(c)
    global pos inStr len
    pos=skip_whitespace(pos,inStr,len);
    if pos > len || inStr(pos) ~= c
        error_pos(sprintf('Expected %c at position %%d', c));
    else
        pos = pos + 1;
        pos=skip_whitespace(pos,inStr,len);
    end
end


function c = next_char
    global pos inStr len
    pos=skip_whitespace(pos,inStr,len);
    if pos > len
        c = [];
    else
        c = inStr(pos);
    end
end


function newpos=skip_whitespace(pos,inStr,len)
    newpos=pos;
    while newpos <= len && isspace(inStr(newpos))
        newpos = newpos + 1;
    end
end


function object = parse_object(varargin)
    parse_char('{');
    object = [];
    if next_char ~= '}'
        while 1
            str = parseStr(varargin{:});
            if isempty(str)
                error_pos('Name of value at position %d cannot be empty');
            end
            parse_char(':');
            val = parse_value(varargin{:});
            object.(valid_field(str))=val;
            if next_char == '}'
                break;
            end
            parse_char(',');
        end
    end
    parse_char('}');
    if(isstruct(object))
        object=struct2jdata(object);
    end
end


function object = parse_array(varargin) % JSON array is written in row-major order
global pos inStr isoct
    parse_char('[');
    object = cell(0, 1);
    dim2=[];
    arraydepth=jsonopt('JSONLAB_ArrayDepth_',1,varargin{:});
    pbar=-1;
    if(isfield(varargin{1},'progressbar_'))
        pbar=varargin{1}.progressbar_;
    end

    if next_char ~= ']'
	if(jsonopt('FastArrayParser',1,varargin{:})>=1 && arraydepth>=jsonopt('FastArrayParser',1,varargin{:}))
            [endpos, e1l, e1r]=matching_bracket(inStr,pos);
            arraystr=['[' inStr(pos:endpos)];
            arraystr=regexprep(arraystr,'"_NaN_"','NaN');
            arraystr=regexprep(arraystr,'"([-+]*)_Inf_"','$1Inf');
            arraystr(arraystr==sprintf('\n'))=[];
            arraystr(arraystr==sprintf('\r'))=[];
            %arraystr=regexprep(arraystr,'\s*,',','); % this is slow,sometimes needed
            if(~isempty(e1l) && ~isempty(e1r)) % the array is in 2D or higher D
        	astr=inStr((e1l+1):(e1r-1));
        	astr=regexprep(astr,'"_NaN_"','NaN');
        	astr=regexprep(astr,'"([-+]*)_Inf_"','$1Inf');
        	astr(astr==sprintf('\n'))=[];
        	astr(astr==sprintf('\r'))=[];
        	astr(astr==' ')='';
        	if(isempty(find(astr=='[', 1))) % array is 2D
                    dim2=length(sscanf(astr,'%f,',[1 inf]));
        	end
            else % array is 1D
        	astr=arraystr(2:end-1);
        	astr(astr==' ')='';
        	[obj, count, errmsg, nextidx]=sscanf(astr,'%f,',[1,inf]);
        	if(nextidx>=length(astr)-1)
                    object=obj;
                    pos=endpos;
                    parse_char(']');
                    return;
        	end
            end
            if(~isempty(dim2))
        	astr=arraystr;
        	astr(astr=='[')='';
        	astr(astr==']')='';
        	astr(astr==' ')='';
        	[obj, count, errmsg, nextidx]=sscanf(astr,'%f,',inf);
        	if(nextidx>=length(astr)-1)
                    object=reshape(obj,dim2,numel(obj)/dim2)';
                    pos=endpos;
                    parse_char(']');
                    if(pbar>0)
                        waitbar(pos/length(inStr),pbar,'loading ...');
                    end
                    return;
        	end
            end
            arraystr=regexprep(arraystr,'\]\s*,','];');
	else
            arraystr='[';
	end
        try
           if(isoct && regexp(arraystr,'"','once'))
                error('Octave eval can produce empty cells for JSON-like input');
           end
           object=eval(arraystr);
           pos=endpos;
        catch
         while 1
            newopt=varargin2struct(varargin{:},'JSONLAB_ArrayDepth_',arraydepth+1);
            val = parse_value(newopt);
            object{end+1} = val;
            if next_char == ']'
                break;
            end
            parse_char(',');
         end
        end
    end
    if(jsonopt('SimplifyCell',0,varargin{:})==1)
      try
        oldobj=object;
        object=cell2mat(object')';
        if(iscell(oldobj) && isstruct(object) && numel(object)>1 && jsonopt('SimplifyCellArray',1,varargin{:})==0)
            object=oldobj;
        elseif(size(object,1)>1 && ismatrix(object))
            object=object';
        end
      catch
      end
    end
    parse_char(']');
    
    if(pbar>0)
        waitbar(pos/length(inStr),pbar,'loading ...');
    end
end


function str = parseStr(varargin)
    global pos inStr len  esc index_esc len_esc
 % len, ns = length(inStr), keyboard
    if inStr(pos) ~= '"'
        error_pos('String starting with " expected at position %d');
    else
        pos = pos + 1;
    end
    str = '';
    while pos <= len
        while index_esc <= len_esc && esc(index_esc) < pos
            index_esc = index_esc + 1;
        end
        if index_esc > len_esc
            str = [str inStr(pos:len)];
            pos = len + 1;
            break;
        else
            str = [str inStr(pos:esc(index_esc)-1)];
            pos = esc(index_esc);
        end
        nstr = length(str);
        switch inStr(pos)
            case '"'
                pos = pos + 1;
                if(~isempty(str))
                    if(strcmp(str,'_Inf_'))
                        str=Inf;
                    elseif(strcmp(str,'-_Inf_'))
                        str=-Inf;
                    elseif(strcmp(str,'_NaN_'))
                        str=NaN;
                    end
                end
                return;
            case '\'
                if pos+1 > len
                    error_pos('End of file reached right after escape character');
                end
                pos = pos + 1;
                switch inStr(pos)
                    case {'"' '\' '/'}
                        str(nstr+1) = inStr(pos);
                        pos = pos + 1;
                    case {'b' 'f' 'n' 'r' 't'}
                        str(nstr+1) = sprintf(['\' inStr(pos)]);
                        pos = pos + 1;
                    case 'u'
                        if pos+4 > len
                            error_pos('End of file reached in escaped unicode character');
                        end
                        str(nstr+(1:6)) = inStr(pos-1:pos+4);
                        pos = pos + 5;
                end
            otherwise % should never happen
                str(nstr+1) = inStr(pos);
                keyboard;
                pos = pos + 1;
        end
    end
    error_pos('End of file while expecting end of inStr');
end
%%-------------------------------------------------------------------------

function num = parse_number(varargin)
    global pos inStr isoct
    currstr=inStr(pos:min(pos+30,end));
    if(isoct~=0)
        numstr=regexp(currstr,'^\s*-?(?:0|[1-9]\d*)(?:\.\d+)?(?:[eE][+\-]?\d+)?','end');
        [num] = sscanf(currstr, '%f', 1);
        delta=numstr+1;
    else
        [num, one, err, delta] = sscanf(currstr, '%f', 1);
        if ~isempty(err)
            error_pos('Error reading number at position %d');
        end
    end
    pos = pos + delta-1;
end
%%-------------------------------------------------------------------------

function val = parse_value(varargin)
    global pos inStr len
    
    if(isfield(varargin{1},'progressbar_'))
        waitbar(pos/len,varargin{1}.progressbar_,'loading ...');
    end
    
    switch(inStr(pos))
        case '"'
            val = parseStr(varargin{:});
            return;
        case '['
            val = parse_array(varargin{:});
            return;
        case '{'
            val = parse_object(varargin{:});
            return;
        case {'-','0','1','2','3','4','5','6','7','8','9'}
            val = parse_number(varargin{:});
            return;
        case 't'
            if pos+3 <= len && strcmpi(inStr(pos:pos+3), 'true')
                val = true;
                pos = pos + 4;
                return;
            end
        case 'f'
            if pos+4 <= len && strcmpi(inStr(pos:pos+4), 'false')
                val = false;
                pos = pos + 5;
                return;
            end
        case 'n'
            if pos+3 <= len && strcmpi(inStr(pos:pos+3), 'null')
                val = [];
                pos = pos + 4;
                return;
            end
    end
    error_pos('Value expected at position %d');
end
%%-------------------------------------------------------------------------

function error_pos(msg)
    global pos inStr len
    poShow = max(min([pos-15 pos-1 pos pos+20],len),1);
    if poShow(3) == poShow(2)
        poShow(3:4) = poShow(2)+[0 -1];  % display nothing after
    end
    msg = [sprintf(msg, pos) ': ' ...
    inStr(poShow(1):poShow(2)) '<error>' inStr(poShow(3):poShow(4)) ];
    error( ['JSONparser:invalidFormat: ' msg] );
end


function str = valid_field(str)
global isoct
% From MATLAB doc: field names must begin with a letter, which may be
% followed by any combination of letters, digits, and underscores.
% Invalid characters will be converted to underscores, and the prefix
% "x0x[Hex code]_" will be added if the first character is not a letter.
    pos=regexp(str,'^[^A-Za-z]','once');
    if(~isempty(pos))
        if(~isoct)
            str=regexprep(str,'^([^A-Za-z])','x0x${sprintf(''%X'',unicode2native($1))}_','once');
        else
            str=sprintf('x0x%X_%s',char(str(1)),str(2:end));
        end
    end
    if(isempty(regexp(str,'[^0-9A-Za-z_]', 'once' )))
        return;
    end
    if(~isoct)
        str=regexprep(str,'([^0-9A-Za-z_])','_0x${sprintf(''%X'',unicode2native($1))}_');
    else
        pos=regexp(str,'[^0-9A-Za-z_]');
        if(isempty(pos))
            return;
        end
        str0=str;
        pos0=[0 pos(:)' length(str)];
        str='';
        for i=1:length(pos)
            str=[str str0(pos0(i)+1:pos(i)-1) sprintf('_0x%X_',str0(pos(i)))];
        end
        if(pos(end)~=length(str))
            str=[str str0(pos0(end-1)+1:pos0(end))];
        end
    end
    %str(~isletter(str) & ~('0' <= str & str <= '9')) = '_';
end


function newdata=struct2jdata(data,varargin)
    %
    % newdata=struct2jdata(data,opt,...)
    %
    % convert a JData object (in the form of a struct array) into an array
    %
    % authors:Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
    %
    % input:
    %      data: a struct array. If data contains JData keywords in the first
    %            level children, these fields are parsed and regrouped into a
    %            data object (arrays, trees, graphs etc) based on JData 
    %            specification. The JData keywords are
    %               "_ArrayType_", "_ArraySize_", "_ArrayData_"
    %               "_ArrayIsSparse_", "_ArrayIsComplex_"
    %      opt: (optional) a list of 'Param',value pairs for additional options 
    %           The supported options include
    %               'Recursive', if set to 1, will apply the conversion to 
    %                            every child; 0 to disable
    %
    % output:
    %      newdata: the covnerted data if the input data does contain a JData 
    %               structure; otherwise, the same as the input.
    %
    % examples:
    %      obj=struct('_ArrayType_','double','_ArraySize_',[2 3],
    %                 '_ArrayIsSparse_',1 ,'_ArrayData_',null);
    %      ubjdata=struct2jdata(obj);
    %
    % license:
    %     BSD License, see LICENSE_BSD.txt files for details 
    %
    % -- this function is part of JSONLab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
    %

    fn=fieldnames(data);
    newdata=data;
    len=length(data);
    if(jsonopt('Recursive',0,varargin{:})==1)
      for i=1:length(fn) % depth-first
        for j=1:len
            if(isstruct(getfield(data(j),fn{i})))
                newdata(j)=setfield(newdata(j),fn{i},jstruct2array(getfield(data(j),fn{i})));
            end
        end
      end
    end
    if(~isempty(strmatch('x0x5F_ArrayType_',fn)) && ~isempty(strmatch('x0x5F_ArrayData_',fn)))
      newdata=cell(len,1);
      for j=1:len
        ndata=cast(data(j).x0x5F_ArrayData_,data(j).x0x5F_ArrayType_);
        iscpx=0;
        if(~isempty(strmatch('x0x5F_ArrayIsComplex_',fn)))
            if(data(j).x0x5F_ArrayIsComplex_)
               iscpx=1;
            end
        end
        if(~isempty(strmatch('x0x5F_ArrayIsSparse_',fn)))
            if(data(j).x0x5F_ArrayIsSparse_)
                if(~isempty(strmatch('x0x5F_ArraySize_',fn)))
                    dim=double(data(j).x0x5F_ArraySize_);
                    if(iscpx && size(ndata,2)==4-any(dim==1))
                        ndata(:,end-1)=complex(ndata(:,end-1),ndata(:,end));
                    end
                    if isempty(ndata)
                        % All-zeros sparse
                        ndata=sparse(dim(1),prod(dim(2:end)));
                    elseif dim(1)==1
                        % Sparse row vector
                        ndata=sparse(1,ndata(:,1),ndata(:,2),dim(1),prod(dim(2:end)));
                    elseif dim(2)==1
                        % Sparse column vector
                        ndata=sparse(ndata(:,1),1,ndata(:,2),dim(1),prod(dim(2:end)));
                    else
                        % Generic sparse array.
                        ndata=sparse(ndata(:,1),ndata(:,2),ndata(:,3),dim(1),prod(dim(2:end)));
                    end
                else
                    if(iscpx && size(ndata,2)==4)
                        ndata(:,3)=complex(ndata(:,3),ndata(:,4));
                    end
                    ndata=sparse(ndata(:,1),ndata(:,2),ndata(:,3));
                end
            end
        elseif(~isempty(strmatch('x0x5F_ArraySize_',fn)))
            if(iscpx && size(ndata,2)==2)
                 ndata=complex(ndata(:,1),ndata(:,2));
            end
            ndata=reshape(ndata(:),data(j).x0x5F_ArraySize_);
        end
        newdata{j}=ndata;
      end
      if(len==1)
          newdata=newdata{1};
      end
    end
end


function endpos = matching_quote(str,pos)
len=length(str);
while(pos<len)
    if(str(pos)=='"')
        if(~(pos>1 && str(pos-1)=='\'))
            endpos=pos;
            return;
        end        
    end
    pos=pos+1;
end
error('unmatched quotation mark');
end

%%-------------------------------------------------------------------------
function [endpos, e1l, e1r, maxlevel] = matching_bracket(str,pos)
    global arraytoken
    level=1;
    maxlevel=level;
    endpos=0;
    bpos=arraytoken(arraytoken>=pos);
    tokens=str(bpos);
    len=length(tokens);
    pos=1;
    e1l=[];
    e1r=[];
    while(pos<=len)
        c=tokens(pos);
        if(c==']')
            level=level-1;
            if(isempty(e1r))
                e1r=bpos(pos);
            end
            if(level==0)
                endpos=bpos(pos);
                return
            end
        end
        if(c=='[')
            if(isempty(e1l))
                e1l=bpos(pos);
            end
            level=level+1;
            maxlevel=max(maxlevel,level);
        end
        if(c=='"')
            pos=matching_quote(tokens,pos+1);
        end
        pos=pos+1;
    end
    if(endpos==0) 
        error('unmatched "]"');
    end
end


function result = GetGoogleSpreadsheet(DOCID)
    % result = GetGoogleSpreadsheet(DOCID)
    % Download a google spreadsheet as csv and import into a Matlab cell array.
    %
    % [DOCID] see the value after 'key=' in your spreadsheet's url
    %           e.g. '0AmQ013fj5234gSXFAWLK1REgwRW02hsd3c'
    %
    % [result] cell array of the the values in the spreadsheet
    %
    % IMPORTANT: The spreadsheet must be shared with the "anyone with the link" option
    %
    % This has no error handling and has not been extensively tested.
    % Please report issues on Matlab FX.
    %
    % DM, Jan 2013
    %
    loginURL = 'https://www.google.com'; 
    csvURL = ['https://docs.google.com/spreadsheet/ccc?key=' DOCID '&output=csv&pref=2'];
    %Step 1: go to google.com to collect some cookies
    cookieManager = java.net.CookieManager([], java.net.CookiePolicy.ACCEPT_ALL);
    java.net.CookieHandler.setDefault(cookieManager);
    handler = sun.net.www.protocol.https.Handler;
    connection = java.net.URL([],loginURL,handler).openConnection();
    connection.getInputStream();
    %Step 2: go to the spreadsheet export url and download the csv
    connection2 = java.net.URL([],csvURL,handler).openConnection();
    result = connection2.getInputStream();
    result = char(readstream(result));
    %Step 3: convert the csv to a cell array
    result = parseCsv(result);
end


function data = parseCsv(data)
% splits data into individual lines
data = textscan(data,'%s','whitespace','\n');
data = data{1};
for ii=1:length(data)
   %for each line, split the string into its comma-delimited units
   %the '%q' format deals with the "quoting" convention appropriately.
   tmp = textscan(data{ii},'%q','delimiter',',');
   data(ii,1:length(tmp{1})) = tmp{1};
end
end


function out = readstream(inStream)
%READSTREAM Read all bytes from stream to uint8
%From: http://stackoverflow.com/a/1323535
import com.mathworks.mlwidgets.io.InterruptibleStreamCopier;
byteStream = java.io.ByteArrayOutputStream();
isc = InterruptibleStreamCopier.getInterruptibleStreamCopier();
isc.copyStream(inStream, byteStream);
inStream.close();
byteStream.close();
out = typecast(byteStream.toByteArray', 'uint8'); 
end


function json=savejson(rootname,obj,varargin)
    %
    % json=savejson(rootname,obj,filename)
    %    or
    % json=savejson(rootname,obj,opt)
    % json=savejson(rootname,obj,'param1',value1,'param2',value2,...)
    %
    % convert a MATLAB object (cell, struct or array) into a JSON (JavaScript
    % Object Notation) string
    %
    % author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
    % created on 2011/09/09
    %
    % $Id$
    %
    % input:
    %      rootname: the name of the root-object, when set to '', the root name
    %        is ignored, however, when opt.ForceRootName is set to 1 (see below),
    %        the MATLAB variable name will be used as the root name.
    %      obj: a MATLAB object (array, cell, cell array, struct, struct array,
    %      class instance).
    %      filename: a string for the file name to save the output JSON data.
    %      opt: a struct for additional options, ignore to use default values.
    %        opt can have the following fields (first in [.|.] is the default)
    %
    %        opt.FileName [''|string]: a file name to save the output JSON data
    %        opt.FloatFormat ['%.10g'|string]: format to show each numeric element
    %                         of a 1D/2D array;
    %        opt.ArrayIndent [1|0]: if 1, output explicit data array with
    %                         precedent indentation; if 0, no indentation
    %        opt.ArrayToStruct[0|1]: when set to 0, savejson outputs 1D/2D
    %                         array in JSON array format; if sets to 1, an
    %                         array will be shown as a struct with fields
    %                         "_ArrayType_", "_ArraySize_" and "_ArrayData_"; for
    %                         sparse arrays, the non-zero elements will be
    %                         saved to _ArrayData_ field in triplet-format i.e.
    %                         (ix,iy,val) and "_ArrayIsSparse_" will be added
    %                         with a value of 1; for a complex array, the 
    %                         _ArrayData_ array will include two columns 
    %                         (4 for sparse) to record the real and imaginary 
    %                         parts, and also "_ArrayIsComplex_":1 is added. 
    %        opt.ParseLogical [0|1]: if this is set to 1, logical array elem
    %                         will use true/false rather than 1/0.
    %        opt.SingletArray [0|1]: if this is set to 1, arrays with a single
    %                         numerical element will be shown without a square
    %                         bracket, unless it is the root object; if 0, square
    %                         brackets are forced for any numerical arrays.
    %        opt.SingletCell  [1|0]: if 1, always enclose a cell with "[]" 
    %                         even it has only one element; if 0, brackets
    %                         are ignored when a cell has only 1 element.
    %        opt.ForceRootName [0|1]: when set to 1 and rootname is empty, savejson
    %                         will use the name of the passed obj variable as the 
    %                         root object name; if obj is an expression and 
    %                         does not have a name, 'root' will be used; if this 
    %                         is set to 0 and rootname is empty, the root level 
    %                         will be merged down to the lower level.
    %        opt.Inf ['"$1_Inf_"'|string]: a customized regular expression pattern
    %                         to represent +/-Inf. The matched pattern is '([-+]*)Inf'
    %                         and $1 represents the sign. For those who want to use
    %                         1e999 to represent Inf, they can set opt.Inf to '$11e999'
    %        opt.NaN ['"_NaN_"'|string]: a customized regular expression pattern
    %                         to represent NaN
    %        opt.JSONP [''|string]: to generate a JSONP output (JSON with padding),
    %                         for example, if opt.JSONP='foo', the JSON data is
    %                         wrapped inside a function call as 'foo(...);'
    %        opt.UnpackHex [1|0]: conver the 0x[hex code] output by loadjson 
    %                         back to the string form
    %        opt.SaveBinary [0|1]: 1 - save the JSON file in binary mode; 0 - text mode.
    %        opt.Compact [0|1]: 1- out compact JSON format (remove all newlines and tabs)
    %
    %        opt can be replaced by a list of ('param',value) pairs. The param 
    %        string is equivallent to a field in opt and is case sensitive.
    % output:
    %      json: a string in the JSON format (see http://json.org)
    %
    % examples:
    %      jsonmesh=struct('MeshNode',[0 0 0;1 0 0;0 1 0;1 1 0;0 0 1;1 0 1;0 1 1;1 1 1],... 
    %               'MeshTetra',[1 2 4 8;1 3 4 8;1 2 6 8;1 5 6 8;1 5 7 8;1 3 7 8],...
    %               'MeshTri',[1 2 4;1 2 6;1 3 4;1 3 7;1 5 6;1 5 7;...
    %                          2 8 4;2 8 6;3 8 4;3 8 7;5 8 6;5 8 7],...
    %               'MeshCreator','FangQ','MeshTitle','T6 Cube',...
    %               'SpecialData',[nan, inf, -inf]);
    %      savejson('jmesh',jsonmesh)
    %      savejson('',jsonmesh,'ArrayIndent',0,'FloatFormat','\t%.5g')
    %
    % license:
    %     BSD or GPL version 3, see LICENSE_{BSD,GPLv3}.txt files for details
    %
    % -- this function is part of JSONLab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
    %

    if(nargin==1)
       varname=inputname(1);
       obj=rootname;
       if(isempty(varname)) 
          varname='root';
       end
       rootname=varname;
    else
       varname=inputname(2);
    end
    if(length(varargin)==1 && ischar(varargin{1}))
       opt=struct('filename',varargin{1});
    else
       opt=varargin2struct(varargin{:});
    end
    opt.IsOctave=exist('OCTAVE_VERSION','builtin');
    if(isfield(opt,'norowbracket'))
        warning('Option ''NoRowBracket'' is depreciated, please use ''SingletArray'' and set its value to not(NoRowBracket)');
        if(~isfield(opt,'singletarray'))
            opt.singletarray=not(opt.norowbracket);
        end
    end
    rootisarray=0;
    rootlevel=1;
    forceroot=jsonopt('ForceRootName',0,opt);
    if((isnumeric(obj) || islogical(obj) || ischar(obj) || isstruct(obj) || ...
            iscell(obj) || isobject(obj)) && isempty(rootname) && forceroot==0)
        rootisarray=1;
        rootlevel=0;
    else
        if(isempty(rootname))
            rootname=varname;
        end
    end
    if((isstruct(obj) || iscell(obj))&& isempty(rootname) && forceroot)
        rootname='root';
    end

    whitespaces=struct('tab',sprintf('\t'),'newline',sprintf('\n'),'sep',sprintf(',\n'));
    if(jsonopt('Compact',0,opt)==1)
        whitespaces=struct('tab','','newline','','sep',',');
    end
    if(~isfield(opt,'whitespaces_'))
        opt.whitespaces_=whitespaces;
    end

    nl=whitespaces.newline;

    json=obj2json(rootname,obj,rootlevel,opt);
    if(rootisarray)
        json=sprintf('%s%s',json,nl);
    else
        json=sprintf('{%s%s%s}\n',nl,json,nl);
    end

    jsonp=jsonopt('JSONP','',opt);
    if(~isempty(jsonp))
        json=sprintf('%s(%s);%s',jsonp,json,nl);
    end

    % save to a file if FileName is set, suggested by Patrick Rapin
    filename=jsonopt('FileName','',opt);
    if(~isempty(filename))
        if(jsonopt('SaveBinary',0,opt)==1)
            fid = fopen(filename, 'wb');
            fwrite(fid,json);
        else
            fid = fopen(filename, 'wt');
            fwrite(fid,json,'char');
        end
        fclose(fid);
    end
end


function txt=obj2json(name,item,level,varargin)
    if(iscell(item))
        txt=cell2json(name,item,level,varargin{:});
    elseif(isstruct(item))
        txt=struct2json(name,item,level,varargin{:});
    elseif(ischar(item))
        txt=str2json(name,item,level,varargin{:});
    elseif(isobject(item)) 
        txt=matlabobject2json(name,item,level,varargin{:});
    else
        txt=mat2json(name,item,level,varargin{:});
    end
end

%%-------------------------------------------------------------------------
function txt=cell2json(name,item,level,varargin)
    txt={};
    if(~iscell(item))
            error('input is not a cell');
    end

    dim=size(item);
    if(ndims(squeeze(item))>2) % for 3D or higher dimensions, flatten to 2D for now
        item=reshape(item,dim(1),numel(item)/dim(1));
        dim=size(item);
    end
    len=numel(item);
    ws=jsonopt('whitespaces_',struct('tab',sprintf('\t'),'newline',sprintf('\n'),'sep',sprintf(',\n')),varargin{:});
    padding0=repmat(ws.tab,1,level);
    padding2=repmat(ws.tab,1,level+1);
    nl=ws.newline;
    bracketlevel=~jsonopt('singletcell',1,varargin{:});
    if(len>bracketlevel)
        if(~isempty(name))
            txt={padding0, '"', checkname(name,varargin{:}),'": [', nl}; name=''; 
        else
            txt={padding0, '[', nl};
        end
    elseif(len==0)
        if(~isempty(name))
            txt={padding0, '"' checkname(name,varargin{:}) '": []'}; name=''; 
        else
            txt={padding0, '[]'};
        end
    end
    for i=1:dim(1)
        if(dim(1)>1)
            txt(end+1:end+3)={padding2,'[',nl};
        end
        for j=1:dim(2)
           txt{end+1}=obj2json(name,item{i,j},level+(dim(1)>1)+(len>bracketlevel),varargin{:});
           if(j<dim(2))
               txt(end+1:end+2)={',' nl};
           end
        end
        if(dim(1)>1)
            txt(end+1:end+3)={nl,padding2,']'};
        end
        if(i<dim(1))
            txt(end+1:end+2)={',' nl};
        end
        %if(j==dim(2)) txt=sprintf('%s%s',txt,sprintf(',%s',nl)); end
    end
    if(len>bracketlevel)
        txt(end+1:end+3)={nl,padding0,']'};
    end
    txt = sprintf('%s',txt{:});
end

%%-------------------------------------------------------------------------
function txt=struct2json(name,item,level,varargin)
    txt={};
    if(~isstruct(item))
        error('input is not a struct');
    end
    dim=size(item);
    if(ndims(squeeze(item))>2) % for 3D or higher dimensions, flatten to 2D for now
        item=reshape(item,dim(1),numel(item)/dim(1));
        dim=size(item);
    end
    len=numel(item);
    forcearray= (len>1 || (jsonopt('SingletArray',0,varargin{:})==1 && level>0));
    ws=struct('tab',sprintf('\t'),'newline',sprintf('\n'));
    ws=jsonopt('whitespaces_',ws,varargin{:});
    padding0=repmat(ws.tab,1,level);
    padding2=repmat(ws.tab,1,level+1);
    padding1=repmat(ws.tab,1,level+(dim(1)>1)+forcearray);
    nl=ws.newline;

    if(isempty(item)) 
        if(~isempty(name)) 
            txt={padding0, '"', checkname(name,varargin{:}),'": []'};
        else
            txt={padding0, '[]'};
        end
        txt = sprintf('%s',txt{:});
        return;
    end
    if(~isempty(name)) 
        if(forcearray)
            txt={padding0, '"', checkname(name,varargin{:}),'": [', nl};
        end
    else
        if(forcearray)
            txt={padding0, '[', nl};
        end
    end
    for j=1:dim(2)
      if(dim(1)>1)
          txt(end+1:end+3)={padding2,'[',nl};
      end
      for i=1:dim(1)
        names = fieldnames(item(i,j));
        if(~isempty(name) && len==1 && ~forcearray)
            txt(end+1:end+5)={padding1, '"', checkname(name,varargin{:}),'": {', nl};
        else
            txt(end+1:end+3)={padding1, '{', nl};
        end
        if(~isempty(names))
          for e=1:length(names)
            txt{end+1}=obj2json(names{e},item(i,j).(names{e}),...
                 level+(dim(1)>1)+1+forcearray,varargin{:});
            if(e<length(names))
                txt{end+1}=',';
            end
            txt{end+1}=nl;
          end
        end
        txt(end+1:end+2)={padding1,'}'};
        if(i<dim(1))
            txt(end+1:end+2)={',' nl};
        end
      end
      if(dim(1)>1)
          txt(end+1:end+3)={nl,padding2,']'};
      end
      if(j<dim(2))
          txt(end+1:end+2)={',' nl};
      end
    end
    if(forcearray)
        txt(end+1:end+3)={nl,padding0,']'};
    end
    txt = sprintf('%s',txt{:});
end

%%-------------------------------------------------------------------------
function txt=str2json(name,item,level,varargin)
    txt={};
    if(~ischar(item))
            error('input is not a string');
    end
    item=reshape(item, max(size(item),[1 0]));
    len=size(item,1);
    ws=struct('tab',sprintf('\t'),'newline',sprintf('\n'),'sep',sprintf(',\n'));
    ws=jsonopt('whitespaces_',ws,varargin{:});
    padding1=repmat(ws.tab,1,level);
    padding0=repmat(ws.tab,1,level+1);
    nl=ws.newline;
    sep=ws.sep;

    if(~isempty(name)) 
        if(len>1)
            txt={padding1, '"', checkname(name,varargin{:}),'": [', nl};
        end
    else
        if(len>1)
            txt={padding1, '[', nl};
        end
    end
    for e=1:len
        val=escapejsonstring(item(e,:));
        if(len==1)
            obj=['"' checkname(name,varargin{:}) '": ' '"',val,'"'];
            if(isempty(name))
                obj=['"',val,'"'];
            end
            txt(end+1:end+2)={padding1, obj};
        else
            txt(end+1:end+4)={padding0,'"',val,'"'};
        end
        if(e==len)
            sep='';
        end
        txt{end+1}=sep;
    end
    if(len>1)
        txt(end+1:end+3)={nl,padding1,']'};
    end
    txt = sprintf('%s',txt{:});
end

%%-------------------------------------------------------------------------
function txt=mat2json(name,item,level,varargin)
    if(~isnumeric(item) && ~islogical(item))
            error('input is not an array');
    end
    ws=struct('tab',sprintf('\t'),'newline',sprintf('\n'),'sep',sprintf(',\n'));
    ws=jsonopt('whitespaces_',ws,varargin{:});
    padding1=repmat(ws.tab,1,level);
    padding0=repmat(ws.tab,1,level+1);
    nl=ws.newline;
    sep=ws.sep;

    if(length(size(item))>2 || issparse(item) || ~isreal(item) || ...
       (isempty(item) && any(size(item))) ||jsonopt('ArrayToStruct',0,varargin{:}))
        if(isempty(name))
            txt=sprintf('%s{%s%s"_ArrayType_": "%s",%s%s"_ArraySize_": %s,%s',...
                  padding1,nl,padding0,class(item),nl,padding0,regexprep(mat2str(size(item)),'\s+',','),nl);
        else
            txt=sprintf('%s"%s": {%s%s"_ArrayType_": "%s",%s%s"_ArraySize_": %s,%s',...
                  padding1,checkname(name,varargin{:}),nl,padding0,class(item),nl,padding0,regexprep(mat2str(size(item)),'\s+',','),nl);
        end
    else
        if(numel(item)==1 && jsonopt('SingletArray',0,varargin{:})==0 && level>0)
            numtxt=regexprep(regexprep(matdata2json(item,level+1,varargin{:}),'^\[',''),']','');
        else
            numtxt=matdata2json(item,level+1,varargin{:});
        end
        if(isempty(name))
            txt=sprintf('%s%s',padding1,numtxt);
        else
            if(numel(item)==1 && jsonopt('SingletArray',0,varargin{:})==0)
                txt=sprintf('%s"%s": %s',padding1,checkname(name,varargin{:}),numtxt);
            else
                txt=sprintf('%s"%s": %s',padding1,checkname(name,varargin{:}),numtxt);
            end
        end
        return;
    end
    dataformat='%s%s%s%s%s';

    if(issparse(item))
        [ix,iy]=find(item);
        data=full(item(find(item)));
        if(~isreal(item))
           data=[real(data(:)),imag(data(:))];
           if(size(item,1)==1)
               % Kludge to have data's 'transposedness' match item's.
               % (Necessary for complex row vector handling below.)
               data=data';
           end
           txt=sprintf(dataformat,txt,padding0,'"_ArrayIsComplex_": ','1', sep);
        end
        txt=sprintf(dataformat,txt,padding0,'"_ArrayIsSparse_": ','1', sep);
        if(size(item,1)==1)
            % Row vector, store only column indices.
            txt=sprintf(dataformat,txt,padding0,'"_ArrayData_": ',...
               matdata2json([iy(:),data'],level+2,varargin{:}), nl);
        elseif(size(item,2)==1)
            % Column vector, store only row indices.
            txt=sprintf(dataformat,txt,padding0,'"_ArrayData_": ',...
               matdata2json([ix,data],level+2,varargin{:}), nl);
        else
            % General case, store row and column indices.
            txt=sprintf(dataformat,txt,padding0,'"_ArrayData_": ',...
               matdata2json([ix,iy,data],level+2,varargin{:}), nl);
        end
    else
        if(isreal(item))
            txt=sprintf(dataformat,txt,padding0,'"_ArrayData_": ',...
                matdata2json(item(:)',level+2,varargin{:}), nl);
        else
            txt=sprintf(dataformat,txt,padding0,'"_ArrayIsComplex_": ','1', sep);
            txt=sprintf(dataformat,txt,padding0,'"_ArrayData_": ',...
                matdata2json([real(item(:)) imag(item(:))],level+2,varargin{:}), nl);
        end
    end
    txt=sprintf('%s%s%s',txt,padding1,'}');
end

%%-------------------------------------------------------------------------
function txt=matlabobject2json(name,item,level,varargin)
    if numel(item) == 0 %empty object
        st = struct();
    else
        % "st = struct(item);" would produce an inmutable warning, because it
        % make the protected and private properties visible. Instead we get the
        % visible properties
        propertynames = properties(item);
        for p = 1:numel(propertynames)
            for o = numel(item):-1:1 % aray of objects
                st(o).(propertynames{p}) = item(o).(propertynames{p});
            end
        end
    end
    txt=struct2json(name,st,level,varargin{:});
end

%%-------------------------------------------------------------------------
function txt=matdata2json(mat,level,varargin)
    ws=struct('tab',sprintf('\t'),'newline',sprintf('\n'),'sep',sprintf(',\n'));
    ws=jsonopt('whitespaces_',ws,varargin{:});
    tab=ws.tab;
    nl=ws.newline;

    if(size(mat,1)==1)
        pre='';
        post='';
        level=level-1;
    else
        pre=sprintf('[%s',nl);
        post=sprintf('%s%s]',nl,repmat(tab,1,level-1));
    end

    if(isempty(mat))
        txt='null';
        return;
    end
    floatformat=jsonopt('FloatFormat','%.10g',varargin{:});
    %if(numel(mat)>1)
        formatstr=['[' repmat([floatformat ','],1,size(mat,2)-1) [floatformat sprintf('],%s',nl)]];
    %else
    %    formatstr=[repmat([floatformat ','],1,size(mat,2)-1) [floatformat sprintf(',\n')]];
    %end

    if(nargin>=2 && size(mat,1)>1 && jsonopt('ArrayIndent',1,varargin{:})==1)
        formatstr=[repmat(tab,1,level) formatstr];
    end

    txt=sprintf(formatstr,mat');
    txt(end-length(nl):end)=[];
    if(islogical(mat) && jsonopt('ParseLogical',0,varargin{:})==1)
       txt=regexprep(txt,'1','true');
       txt=regexprep(txt,'0','false');
    end
    %txt=regexprep(mat2str(mat),'\s+',',');
    %txt=regexprep(txt,';',sprintf('],\n['));
    % if(nargin>=2 && size(mat,1)>1)
    %     txt=regexprep(txt,'\[',[repmat(sprintf('\t'),1,level) '[']);
    % end
    txt=[pre txt post];
    if(any(isinf(mat(:))))
        txt=regexprep(txt,'([-+]*)Inf',jsonopt('Inf','"$1_Inf_"',varargin{:}));
    end
    if(any(isnan(mat(:))))
        txt=regexprep(txt,'NaN',jsonopt('NaN','"_NaN_"',varargin{:}));
    end
end

%%-------------------------------------------------------------------------
function newname=checkname(name,varargin)
    isunpack=jsonopt('UnpackHex',1,varargin{:});
    newname=name;
    if(isempty(regexp(name,'0x([0-9a-fA-F]+)_','once')))
        return
    end
    if(isunpack)
        isoct=jsonopt('IsOctave',0,varargin{:});
        if(~isoct)
            newname=regexprep(name,'(^x|_){1}0x([0-9a-fA-F]+)_','${native2unicode(hex2dec($2))}');
        else
            pos=regexp(name,'(^x|_){1}0x([0-9a-fA-F]+)_','start');
            pend=regexp(name,'(^x|_){1}0x([0-9a-fA-F]+)_','end');
            if(isempty(pos))
                return;
            end
            str0=name;
            pos0=[0 pend(:)' length(name)];
            newname='';
            for i=1:length(pos)
                newname=[newname str0(pos0(i)+1:pos(i)-1) char(hex2dec(str0(pos(i)+3:pend(i)-1)))];
            end
            if(pos(end)~=length(name))
                newname=[newname str0(pos0(end-1)+1:pos0(end))];
            end
        end
    end
end

%%-------------------------------------------------------------------------
function newstr=escapejsonstring(str)
    newstr=str;
    isoct=exist('OCTAVE_VERSION','builtin');
    if(isoct)
       vv=sscanf(OCTAVE_VERSION,'%f');
       if(vv(1)>=3.8)
           isoct=0;
       end
    end
    if(isoct)
      escapechars={'\\','\"','\/','\a','\f','\n','\r','\t','\v'};
      for i=1:length(escapechars);
        newstr=regexprep(newstr,escapechars{i},escapechars{i});
      end
      newstr=regexprep(newstr,'\\\\(u[0-9a-fA-F]{4}[^0-9a-fA-F]*)','\$1');
    else
      escapechars={'\\','\"','\/','\a','\b','\f','\n','\r','\t','\v'};
      for i=1:length(escapechars);
        newstr=regexprep(newstr,escapechars{i},regexprep(escapechars{i},'\\','\\\\'));
      end
      newstr=regexprep(newstr,'\\\\(u[0-9a-fA-F]{4}[^0-9a-fA-F]*)','\\$1');
    end
end

%%-------------------------------------------------------------------------
function str = BuildRange(sheetpos, ca)
% for cell array, ca, build the range in A1 notation
% https://developers.google.com/sheets/reference/rest/v4/spreadsheets.values#ValueRange
%
% don't know who to credit the following code to
%
% convert number, number format to alpha, number format
%t = [floor(c/27) + 64 floor((c - 1)/26) - 2 + rem(c - 1, 26) + 65];

% The result can have up to three alpha digits. By Peter
% https://www.mathworks.com/matlabcentral/newsreader/view_thread/85589

r = sheetpos(1);
c = sheetpos(2);
Digits = zeros(1, 3);
% Convert number-number format to alpha-number format.
Digits(1) = max(floor(((c - 1) / 26 - 1) / 26), 0);
Digits(2) = floor((c - Digits(1) * 26 * 26 - 1) / 26);
Digits(3) = rem(c - 1, 26) + 1;
% Delete negative numbers and convert blank cells to spaces.
Digits(Digits > 0) = Digits(Digits > 0) + 64;
Digits(Digits == 0) = Digits(Digits == 0) + 32;
% There may be leading spaces, so trim them away.
start = strtrim([char(Digits), num2str(r)]);

r = sheetpos(1)+size(ca,1)-1;
c = sheetpos(2)+size(ca,2)-1;
Digits = zeros(1, 3);
% Convert number-number format to alpha-number format.
Digits(1) = max(floor(((c - 1) / 26 - 1) / 26), 0);
Digits(2) = floor((c - Digits(1) * 26 * 26 - 1) / 26);
Digits(3) = rem(c - 1, 26) + 1;
% Delete negative numbers and convert blank cells to spaces.
Digits(Digits > 0) = Digits(Digits > 0) + 64;
Digits(Digits == 0) = Digits(Digits == 0) + 32;
% There may be leading spaces, so trim them away.
stop = strtrim([char(Digits), num2str(r)]);

str = [start ':' stop];
end