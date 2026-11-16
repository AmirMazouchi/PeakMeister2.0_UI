# python3 tests/investigate-browser.py /path/to/tiny.pwiz.1.1.mz5 --url http://127.0.0.1:3840
# Requires Playwright and Chromium; run against the actual app.
from playwright.sync_api import sync_playwright
import json
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("fixture",help="ProteoWizard tiny.pwiz.1.1.mz5")
parser.add_argument("--url",default="http://127.0.0.1:3840")
args=parser.parse_args()
with sync_playwright() as p:
    browser=p.chromium.launch(executable_path='/usr/bin/chromium',headless=True,args=['--no-sandbox','--disable-dev-shm-usage'])
    page=browser.new_page(viewport={'width':1600,'height':1100})
    errors=[]
    page.on('pageerror',lambda e:errors.append(str(e)))
    print('Opening',flush=True)
    page.goto(args.url,wait_until='networkidle',timeout=120000)
    page.wait_for_function('window.Shiny && Shiny.shinyapp && Shiny.shinyapp.$socket.readyState === 1')
    assert page.locator('#autorun').is_checked()
    page.locator('a[data-value="engine"]').click()
    page.locator('#autorun').uncheck()
    page.locator('#mz5Files').set_input_files(args.fixture)
    page.locator('a[data-value="investigate"]').click()
    page.wait_for_function('document.querySelector("#investigate_status").textContent.includes("Engine data:")',timeout=120000)
    page.wait_for_function('document.getElementById("investigate_median_plot")?.data?.length > 0',timeout=120000)
    print('Plot ready',flush=True)
    assert page.locator('#investigate_time_min').input_value()=='1.2'
    assert page.locator('#investigate_median_statistic').input_value()=='max'
    assert page.locator('#investigate_raw_file').count()==0
    assert page.evaluate('$("#investigate_window_position").data("ionRangeSlider").metaboLive === true')
    page.evaluate('''() => { var g=document.getElementById('investigate_median_plot');
      window.metaboGraph=g;window.metaboRestyles=[];
      g.on('plotly_restyle',function(e){window.metaboRestyles.push(e)});
      return Plotly.relayout(g,{'xaxis.range':[162,165],'yaxis.range':[-10,200]}); }''')
    page.wait_for_timeout(500)
    # Move from a populated final partial slice to the initial slice (empty in
    # this tiny mz5). Both must retain the original DOM widget and both axes.
    page.locator('#investigate_window_width').evaluate("el => el.selectize.setValue('1')")
    page.wait_for_timeout(500)
    # Capture actual handle geometry from the live ionRangeSlider container.
    page.locator('#investigate_window_position').evaluate("el => $(el).data('ionRangeSlider').$cache.cont[0].scrollIntoView({block:'center'})")
    before_drag=page.evaluate('window.metaboRestyles.length')
    box=page.evaluate('''() => {var s=$('#investigate_window_position').data('ionRangeSlider');
      var r=s.$cache.cont[0].getBoundingClientRect(),h=s.$cache.s_single[0].getBoundingClientRect();return {x:r.x,y:r.y,w:r.width,h:r.height,hx:h.x+h.width/2,hy:h.y+h.height/2};}''')
    page.mouse.move(box['hx'],box['hy'])
    page.mouse.down()
    page.mouse.move(box['x']+box['w']*.8,box['hy'],steps=15)
    page.wait_for_timeout(700)
    during=page.evaluate('({live:Shiny.shinyapp.$inputValues.investigate_window_live,updates:window.metaboRestyles.length})')
    page.mouse.up()
    page.wait_for_timeout(700)
    result=page.evaluate('''() => {var g=document.getElementById('investigate_median_plot');return {
      same:g===window.metaboGraph,x:g.layout.xaxis.range,y:g.layout.yaxis.range,
      interval:document.querySelector('#investigate_window_interval').textContent,
      updates:window.metaboRestyles.length,live:Shiny.shinyapp.$inputValues.investigate_window_live};}''')
    assert result['same'],result
    assert result['x']==[162,165] and result['y']==[-10,200],result
    assert during['updates']>before_drag and during['live'] is not None,(during,result)
    assert not errors,errors
    print(json.dumps({'PASS':'Engine upload and live browser drag; unchanged SVG widget and x/y limits','during_drag':during,'result':result}))
    browser.close()
