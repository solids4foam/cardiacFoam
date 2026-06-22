import os
import re
from pathlib import Path

def find_info_only_params(src_dir):
    lookup_regex = re.compile(r'([A-Za-z0-9_]+)(?:\s*=\s*|\s*\(\s*(?:[A-Za-z0-9_]+\.)?lookup(?:OrDefault)?(?:<[^>]+>)?\s*\(\s*"([^"]+)")')
    
    # We will find all variables that are assigned from a lookup.
    for root, _, files in os.walk(src_dir):
        for f in files:
            if f.endswith('.C'):
                path = os.path.join(root, f)
                with open(path, 'r', encoding='utf-8') as file:
                    content = file.read()
                
                # Simple strip comments
                content_no_comments = re.sub(r'//.*', '', content)
                content_no_comments = re.sub(r'/\*.*?\*/', '', content_no_comments, flags=re.DOTALL)
                
                # find variable initializations from dict lookups
                # this regex matches: varName(dict.lookup("key")) or varName = dict.lookup("key")
                # we also want to catch: type varName(dict.lookup("key"))
                matches = re.finditer(r'([A-Za-z0-9_]+)\s*(?:=|\()\s*(?:[A-Za-z0-9_\.]*)lookup(?:OrDefault)?(?:<[^>]+>)?\s*\(\s*"([^"]+)"', content_no_comments)
                for m in matches:
                    var_name = m.group(1)
                    param_name = m.group(2)
                    
                    # exclude standard types being mistakenly matched
                    if var_name in ['word', 'scalar', 'label', 'bool', 'vector', 'tensor', 'string', 'readScalar', 'readLabel', 'readBool']:
                        continue
                        
                    # Now count occurrences of var_name in the file
                    # We tokenize the file (without comments) to find exact word matches
                    words = re.findall(r'[A-Za-z0-9_]+', content_no_comments)
                    var_count = words.count(var_name)
                    
                    # If it's used very few times, maybe it's only in Info.
                    # Let's check where it's used. We can split the content into statements (by semicolon)
                    statements = content_no_comments.split(';')
                    info_only = True
                    for stmt in statements:
                        # If the statement contains the variable
                        if re.search(r'\b' + re.escape(var_name) + r'\b', stmt):
                            # If it's not the initialization statement (which contains lookup)
                            if 'lookup' not in stmt:
                                # Is it an Info statement?
                                if 'Info' not in stmt and 'Warning' not in stmt:
                                    info_only = False
                                    break
                    
                    if info_only and var_count > 1: # var_count > 1 means it's used at least once after initialization
                        print(f"File: {path}")
                        print(f"  Param: '{param_name}' -> Var: {var_name}")
                        print(f"  Occurrences: {var_count}")

find_info_only_params('src')
