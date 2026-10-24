from b1_common import *
if __name__=='__main__':
    stage=sys.argv[1]
    try:
        if stage=='initialize':initialize()
        elif stage=='diagnostics':
            from diagnostics import run
            run()
        elif stage=='train':
            from train_family import run
            run()
        elif stage=='report':
            from deliver import run
            run()
        elif stage in ['pre_tests','final_tests']:
            import subprocess,xml.etree.ElementTree as ET
            xml=E/'logs'/('pre_model_tests.xml' if stage=='pre_tests' else 'final_tests.xml')
            target=E/'tests/test_rules.py' if stage=='pre_tests' else E/'tests'
            result=subprocess.run([sys.executable,'-m','pytest',str(target),'-q','--junitxml='+str(xml)])
            if result.returncode:raise RuntimeError('Required tests failed')
            root=ET.parse(xml).getroot();s=root if root.tag=='testsuite' else root.find('testsuite')
            if stage=='pre_tests':dump(dict(passed=True,tests=int(s.attrib['tests']),source=str(xml),before_research_fitting=True),E/'logs/pre_model_test_gate.json')
        else:raise ValueError(stage)
    except Exception as e:
        log('stage_failure',stage=stage,error=repr(e),traceback=traceback.format_exc());raise
