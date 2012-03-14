import sqlite3

from flask import Flask, request, redirect, url_for, g, render_template

import compute

db  = compute.filenames

app = Flask(__name__)


def query_html(s="N<=100 and k=2 and prec>0"):
    return \
"""<form action="spaces" method="get">
  SQL query on N, k, i, newforms, prec:<br><input type="text" name="query" size=120 value="%s"/> 
  <input type="submit" value="Submit" />
</form>"""%s

def space_html(N, k, i):
    return "<i>S</i><sub>%s</sub>(%s%s)"%(k, N,'' if i==0 else ', &chi;<sub>%s</sub>'%i)

def newform_html(N,k,i,j, in_space=False):
    return '<i>f</i><sub>%s</sub>'%j + ('&isin; %s'%space_html(N,k,i) if in_space else '')

def space_url(N,k,i):
    return '/space?N=%s&k=%s&i=%s'%(N,k,i)

def newform_url(N,k,i,j):
    return '/newform?N=%s&k=%s&i=%s&j=%s'%(N,k,i,j)

@app.route("/")
def main_page():
    return query_html()

@app.route("/spaces", methods=['GET'])
def spaces():
    query = request.args.get('query','')
    if ';' in query: query = query.split(';')[0]

    v = []
    max_lines = 1000
    try:
        result = db.known(query)
    except sqlite3.OperationalError, msg:
        return query_html(query) + '<br>%s'%msg
    num = 0
    result = list(result)
    num = len(result)
    for (N, k, i, newforms, maxp) in result:
        line = '<a style="text-decoration: none;" href="%s">%s</a>'%(
            space_url(N,k,i), space_html(N,k,i))
        line += ' &nbsp;--&nbsp; '        
        if newforms == -1:
            line += '(newforms not computed)'
        else:
            for j in range(newforms):
                line += '<a style="text-decoration: none;" href="%s">%s</a>&nbsp;&nbsp;&nbsp;'%(
                    newform_url(N,k,i,j), newform_html(N,k,i,j))
        v.append(line)
        if len(v) >= max_lines:
            v.insert(0, "<b>Showing first %s results</b><br>"%len(v))
            break
    return query_html(query) + '<h1>Spaces (%s listed below)</h1>'%num + '<br>'.join(v)

@app.route("/space", methods=['GET'])
def space():
    N = int(request.args.get('N', 37))
    k = int(request.args.get('k', 2))
    i = int(request.args.get('i', 0))
    query = 'N=%s and k=%s and i=%s'%(N,k,i)
    v = list(db.known(query))
    if len(v) == 0:
        return "No information in database about %s"%space_html(N,k,i)
    N, k, i, newforms, maxp = v[0]
    if newforms == -1:
        return "Newform data not yet computed in %s"%space_html(N,k,i)
    if newforms == 0:
        return "No newforms in %s"%space_html(N,k,i)
    w = []
    for j in range(newforms):
        line = '<a style="text-decoration: none;" href="%s">%s</a>'%(
            newform_url(N,k,i,j), newform_html(N,k,i,j))
        w.append(line)
    return query_html(query) + '<h1>Newforms</h1>' + '<br>'.join(w)

@app.route("/newform", methods=['GET'])
def newform():
    N = int(request.args.get('N', 37))
    k = int(request.args.get('k', 2))
    i = int(request.args.get('i', 0))
    j = int(request.args.get('j', 0))
    s = query_html('N=%s and k=%s and i=%s'%(N,k,i))
    s += '<h1>%s</h1>'%newform_html(N,k,i,j, in_space=True)

    # this could take significant time
    A = compute.load_factor(N,k,i,j)
    s += 'Newform of degree %s<br>'%A.dimension()

    s += '<a href="/aplist_powerbasis?N=%s&k=%s&i=%s&j=%s">power basis</a>'%(
        N,k,i,j)
    
    return s


@app.route("/aplist_powerbasis", methods=['GET'])
def aplist_powerbasis():
    N = int(request.args.get('N', 37))
    k = int(request.args.get('k', 2))
    i = int(request.args.get('i', 0))
    j = int(request.args.get('j', 0))
    v = compute.load(db.factor_dual_eigenvector(N,k,i,j,makedir=False))
    aplist = compute.load(db.factor_aplist(N,k,i,j,False,100))
    K = v.parent().base_ring()
    s = '<pre>\n'
    if hasattr(K, 'defining_polynomial'):
        s += 'f = %s;\n'%str(K.defining_polynomial()).replace(' ', '')
    ap = str(list(aplist*v)).replace(' ','')
    if K.absolute_degree() >= 4:
        ap = ap.replace(',',',\n')
    s += 'ap = %s'%(ap)
    s += '\n</pre>'
    return s

if __name__ == '__main__':
    app.run(debug=True, port=5389)
