def file_work(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    
    print(lines[0])
    lines_out = []
    for line in lines:
        line_out = ''
        for i in range(len(line)-1):
            if line[i] != ' ':
                line_out += line[i]
            elif line[i] == ' ' and line[i-1] != ' ' and line[i-1] != ',':
                line_out += ','
            else:
                pass
        lines_out.append(line_out + '\n')
    
    print(lines_out[0])
    file_out = file[:-4] + '_out.txt'
    with open(file_out, 'w') as out:
        out.writelines(lines_out)
    
    print('Done!')

if __name__ == '__main__':
    file_work('EnvironmentMeasuresv06.txt')
        

