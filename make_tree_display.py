
import pickle
import CasSites
import os
import globals


def sub_tree_display(path, results_url, candidates_lst, f, consider_homology, counter, genes_names, genes_lst,
                     repetitions_flag, genomeAssembly, use_thr):
    """

    :param path:
    :param results_url:
    :param candidates_lst:
    :param f:
    :param consider_homology:
    :param counter:
    :param genes_names:
    :param genes_lst:
    :param repetitions_flag:
    :param genomeAssembly:
    :param use_thr:
    """
    runId = ''.join(list(filter(str.isdigit, path)))

    if consider_homology:
        f.write("<table id=" + str(counter) + " class=\"table table-striped table-hover\" style='display:none; width:100%'>\n")
    else:
        f.write("<br><u><h2>")

        if use_thr == 0:
            f.write("A single sgRNA that could best target the entire gene set")
        else:
            if not results_url:
                f.write("The best set of sgRNAs that could target the genes set with scores above the input threshold")
            else:
                f.write("A single sgRNA optimized to target most of the genes, each gene with a score above the input threshold")

        f.write("</h2></u>\n")
        f.write("<p><u><b><font style=\"font-size:18px;\">Filters:</font></b></u></p>")
        f.write(
            "<font style=\"font-size:16px;\">&#9658; Show only targets with up to <select id=\"selectMM\" onchange=\"javascript:filterMM();\">\n")
        f.write("\t<option value=\"" + str(100) + "\">""</option>\n")  # new
        for i in range(1, 21):  # new
            f.write("\t<option value=\"" + str(i) + "\">" + str(i) + "</option>\n")
        f.write("</select>\n mismatches\n")

        if results_url:
            if repetitions_flag:
                repetitions_str = "For each subset of targets, multiple sgRNAs are designed, but <u>only the best one</u> is presented in the table below. To see all the sgRNAs"
            else:
                repetitions_str = "For each subset of targets, multiple sgRNAs are designed. To see </u>only the best sgRNA</u> for each subset"

            f.write(
                "<br><br>&#9658; " + repetitions_str + " <a id=\"abc\" href=\"" + results_url + "\" style=\"font-size:16px; text-decoration: underline;\">click here</a>.\n")

        f.write("<input type=\"hidden\" id=\"scrl\" name=\"user\" value=\"0\" />\n")

        if not consider_homology and not results_url:
            crista_url = "http://crista.tau.ac.il/results/" + runId + "_crispys_greedy/output.php"
        else:
            crista_url = "http://crista.tau.ac.il/results/" + runId + "_crispys/output.php"

        if genomeAssembly != "Not_selected":
            f.write(
                "<br><br>CRISTA*: For the genome wide-wide off-targets detection results for the best sgRNAs click <a href=\"" + crista_url + "\" target=\"_blank\" style=\"font-size:16px; text-decoration: underline;\">here</a>.</font><br><br>\n")
        else:
            f.write("</font><br><br><br>\n")

        f.write("<table class=\"table\">\n")

    header_row = "<tr>\n\t<th></th>\n\t<th>sgRNA</th>\n\t<th>Score</th>\n\t<th>Genes</th>\n\t<th>Genes score</th>\n\t<th>Target site</th>\n\t<th>#mms</th>\n\t<th>Position</th>\n\t<th>Find off-targets</th>\n</tr>\n"  # new
    f.write(header_row)

    n_candidate = 0
    if consider_homology:
        n_use_crista_link = 1
    else:
        n_use_crista_link = 5

    for candidate in candidates_lst:
        n_candidate += 1
        if candidate.off_targets and genomeAssembly != "Not_selected":
            crista_label = "CRISTA*"
            if not consider_homology and not results_url:
                crista_url = "http://crista.tau.ac.il/results/" + runId + "_crispys_greedy/output.php"
            else:
                crista_url = "http://crista.tau.ac.il/results/" + runId + "_crispys/output.php"
        else:
            crista_label = "CRISTA"
            crista_url = "http://crista.tau.ac.il/findofftargets.html#form/sgRNA_seq=" + candidate.seq

        f.write(
            "<tr class='top_row'>\n\t<td></td>\n\t<td></td>\n\t<td></td>\n\t<td></td>\n\t<td></td>\n\t<td></td>\n\t<td></td>\n\t<td></td>\n\t<td></td>\n</tr>\n")
        f.write("<tr>\n")
        num_of_targets = 0
        for targets in candidate.targets_dict.values():
            num_of_targets += len(targets)
        first_gene = 1
        l = list(candidate.targets_dict.items())
        l.sort(key=lambda item: candidate.genes_score_dict[item[0]], reverse=True)

        for gene, targets in l:
            targets.sort(key=lambda target: len(target[1]))

            gene_seq = genes_lst[genes_names.index(gene)]
            seen_sites = dict()
            first_target = 1
            for target in targets:
                pos = find_pos(target, gene_seq, seen_sites)
                if first_target == 1 and first_gene == 1:
                    f.write(
                        "\t<td class='" + candidate.seq + "' style='display:none' onclick='hide(\"" + candidate.seq + "\")'>-</td>\n")
                    f.write("\t<td class='" + candidate.seq + "' onclick='show_row(\"" + candidate.seq + "\")'>+</td>\n")
                    f.write("\t<td class='" + candidate.seq + "' style='display:none'>" + candidate.seq + "</td>\n")
                    f.write("\t<td class='" + candidate.seq + "'>" + candidate.seq + "</td>\n")
                    f.write(
                        "\t<td class='" + candidate.seq + "' style='display:none'>" + str(candidate.cut_expectation)[:5] + "</td>\n")
                    f.write("\t<td class='" + candidate.seq + "'>" + str(candidate.cut_expectation)[:5] + "</td>\n")
                    f.write("\t<td class='" + candidate.seq + "' style='display:none'>" + gene + "</td>\n")
                    f.write("\t<td class='" + candidate.seq + "'>" + gene + "</td>\n")
                    score = str(candidate.genes_score_dict[gene])
                    if len(score) > 5:
                        score = score[:5]
                    f.write("\t<td class='" + candidate.seq + "' style='display:none'>" + score + "</td>\n")
                    f.write("\t<td class='" + candidate.seq + "'>" + score + "</td>\n")
                    f.write("\t<td style='font-family:Courier new'>" + change_mismatch_to_lowercase(target[0], target[1].keys()) + "</td>\n")
                    f.write("\t<td id=\"mms\">" + str(len(target[1])) + "</td>\n")
                    f.write("\t<td>" + pos + "</td>\n")
                    f.write(
                        "\t<td class='" + candidate.seq + "' style='display:none'><a target=\"_blank\" href=\"http://tefor.net/crispor/crispor.cgi?seq=" + candidate.seq + "&pam=NGG\">CRISPOR</a> | <a target=\"_blank\" href=\"" + crista_url + "\">" + crista_label + "</a></td>\n")
                    f.write(
                        "\t<td class='" + candidate.seq + "'><a target=\"_blank\" href=\"http://tefor.net/crispor/crispor.cgi?&seq=" + candidate.seq + "&pam=NGG\">CRISPOR</a> | <a target=\"_blank\" href=\"" + crista_url + "\">" + crista_label + "</a></td>\n")
                    f.write("</tr>\n")
                    first_target = 0
                    continue
                if first_target != 1:
                    f.write("<tr class='" + candidate.seq + "_r' style='display:none'>\n")
                    f.write("\t<td></td>\n")
                    f.write("\t<td></td>\n")
                    f.write("\t<td></td>\n")
                    f.write("\t<td></td>\n")
                    f.write("\t<td></td>\n")
                    f.write("\t<td style='font-family:Courier new'>" + change_mismatch_to_lowercase(target[0], target[1].keys()) + "</td>\n")
                    f.write("\t<td id=\"mms\">" + str(len(target[1])) + "</td>\n")  # new
                    f.write("\t<td>" + pos + "</td>\n")  # new
                    f.write("\t<td></td>\n")
                    f.write("</tr>\n")
                if first_target == 1 and first_gene != 1:
                    score = str(candidate.genes_score_dict[gene])
                    if len(score) > 5:
                        score = score[:5]
                    f.write("<tr class='" + candidate.seq + "_r' style='display:none'>\n")
                    f.write("\t<td></td>\n")
                    f.write("\t<td></td>\n")
                    f.write("\t<td></td>\n")
                    f.write("\t<td>" + gene + "</td>\n")
                    f.write("\t<td>" + score + "</td>\n")
                    f.write("\t<td style='font-family:Courier new'>" + change_mismatch_to_lowercase(target[0], target[1].keys()) + "</td>\n")
                    f.write("\t<td id=\"mms\">" + str(len(target[1])) + "</td>\n")  # new
                    f.write("\t<td>" + pos + "</td>\n")  # new
                    f.write("\t<td></td>\n")
                    f.write("</tr>\n")
                    first_target = 0
            first_gene = 0

    f.write(
        "<tr class='top_row'>\n\t<td></td>\n\t<td></td>\n\t<td></td>\n\t<td></td>\n\t<td></td>\n\t<td></td>\n\t<td></td>\n\t<td></td>\n\t<td></td>\n</tr>")
    f.write("</table>\n")


def change_mismatch_to_lowercase(target_str, mm_lst):
    """

    :param target_str:
    :param mm_lst:
    :return:
    """
    target_in_lst = list(target_str)
    for place in mm_lst:
        target_in_lst[place] = target_in_lst[place].lower()
    return ''.join(target_in_lst)


def find_pos(target, gene_sequence, seen_sites):
    """sgRNA_targets is a list of target sites
	returns"""

    target_seq = target[0]
    if target_seq in seen_sites:
        directions_lst = seen_sites[target_seq]
    else:
        directions_lst = [0, 0]
    position = gene_sequence.find(target_seq, directions_lst[0])
    if position != -1:
        update_seen_sites_dict(seen_sites, target_seq, 0, position)
        position = str(position) + '+'
    else:
        position = CasSites.give_complementary(gene_sequence).find(target_seq, directions_lst[1])
        update_seen_sites_dict(seen_sites, target_seq, 1, position)
        position = str(position) + '-'
    if position == -1:
        position = ''
    return position


def update_seen_sites_dict(d, site_seq, direction, position):
    """
	d: dict: key: site_seq; val: [num of occurrences, directions_lst]
	direction: 0 for sense, 1 for antisense
	"""
    if site_seq in d:
        directions_lst = d[site_seq]
        directions_lst[direction] = position + 1
        d[site_seq] = directions_lst
    else:
        directions_lst = [0, 0]
        directions_lst[direction] = position + 1
        d[site_seq] = directions_lst


def tree_display(path, consider_homology=False, data='regular', genomeAssembly='Not_selected', use_thr=0):
    """

    :param path:
    :param consider_homology: :
    :param data:
    :param genomeAssembly:
    :param use_thr:
    """
    multicrispr_url = "http://multicrispr.tau.ac.il"
    runId = ''.join(list(filter(str.isdigit, path)))
    results_url = ""

    if data == 'greedy_set_cover':
        candidates_lst = pickle.load(open(path + "/greedy_cover.p", "rb"))
        filepath = path + "/cover_table.html"
        f = open(filepath, 'w')
        use_thr = True
    elif data == 'removed_repetitions':
        candidates_lst = pickle.load(open(path + '/res_in_lst_removed_rep.p', "rb"))
        filepath = path + "/removed_repetitions_table.html"
        f = open(filepath, 'w')
        results_url = multicrispr_url + "/results.html?jobId=" + runId + "&removeRepetitions=0"
    else:
        candidates_lst = pickle.load(open(path + "/res_in_lst.p", "rb"))
        filepath = path + "/the_table.html"
        f = open(filepath, 'w')
        results_url = multicrispr_url + "/results.html?jobId=" + runId + "&removeRepetitions=1"
    genes_names = pickle.load(open(path + "/genes_names_list.p", "rb"))
    genes_list = pickle.load(open(path + '/genes_list.p', 'rb'))

    # create crispys page
    if genomeAssembly != "Not_selected":

        n_candidate = 0
        sgrna_list = ''

        # create sgrna_list
        if not consider_homology:
            for c in candidates_lst:
                if c.off_targets:
                    if sgrna_list:
                        sgrna_list += ','
                    sgrna_list += c.seq
        else:
            for subgroup_item in candidates_lst:
                for c in subgroup_item.candidate_list:
                    if c.off_targets:
                        if sgrna_list:
                            sgrna_list += ','
                        sgrna_list += c.seq

        # write and execute shell command
        if sgrna_list:
            cmd = "python /bioseq/crista/CRISTA_online/multiple_sgRNAs_offtargets_search.py -s " + sgrna_list + " -g " + genomeAssembly + " -m 1 -n " + runId + "_crispys"
            if data == 'greedy_set_cover' and not consider_homology:
                cmd += "_greedy"
            os.system(f"{globals.ssh_connect} {cmd}\"")

    counter = 0

    if not consider_homology:
        sub_tree_display(path, results_url, candidates_lst, f, consider_homology, counter, genes_names, genes_list,
                         data == 'removed_repetitions', genomeAssembly, use_thr)
    else:
        f.write(
            "<font style=\"font-size:18px;\"><b>The designed sgRNAs for the genes in your input are listed in the table below.<br>Every section of the table corresponds to a subgroup of homologous genes as specified by the internal nodes of the constructed genes tree.<br>The name of the subgroup and the list of genes are given in the header of each section.</b></font><br><br>")

        f.write("<br><u><h2>")

        if use_thr == 0:
            f.write("Multiple sgRNAs, each optimized towards subgroup of homologous genes")
        else:
            f.write(
                "Multiple sgRNAs, each optimized towards subgroup of homologous genes with a score above the input threshold")

        f.write("</h2></u>\n")
        f.write("<p><u><b><font style=\"font-size:18px;\">Filters:</font></b></u></p>")
        f.write(
            "<font style=\"font-size:16px;\">&#9658; Show only targets with up to <select id=\"selectMM\" onchange=\"javascript:filterMM();\">\n")  # new

        f.write("\t<option value=\"" + str(100) + "\">""</option>\n")  # new
        for i in range(1, 21):  # new
            f.write("\t<option value=\"" + str(i) + "\">" + str(i) + "</option>\n")  # new
        f.write("</select>\n mismatches\n")  # new

        if results_url:
            if data == 'removed_repetitions':
                repetitions_str = "For each subset of targets, multiple sgRNAs are designed, but <u>only the best one</u> is presented in the table below. To see all the sgRNAs"
            else:
                repetitions_str = "For each subset of targets, multiple sgRNAs are designed. To see </u>only the best sgRNA</u> for each subset"

            f.write(
                "<br><br>&#9658; " + repetitions_str + " <a id=\"abc\" href=\"" + results_url + "\" style=\"font-size:16px; text-decoration: underline;\">click here</a>.\n")

        f.write("<input type=\"hidden\" id=\"scrl\" name=\"user\" value=\"0\" />\n")

        crista_url = "http://crista.tau.ac.il/results/" + runId + "_crispys/output.php"

        if genomeAssembly != "Not_selected":
            f.write(
                "<br><br>CRISTA*: For the genome wide-wide off-targets detection results for the best sgRNAs click <a href=\"" + crista_url + "\" target=\"_blank\" style=\"font-size:16px; text-decoration: underline;\">here</a>.</font><br><br>\n")
        else:
            f.write("</font><br><br><br>\n")

        for subgroup_item in candidates_lst:
            # create the main table
            counter += 1
            genes = str(subgroup_item.genes_lst)[1:-1].replace("'", "")

            f.write("<table class=\"table group\" style=\"font-size:20px; font-family:aharoni; \">\n")
            f.write("  <tr class=\"group\">\n       <td onclick='show(this, " + str(
                counter) + ")'>+ " + subgroup_item.name + ": " + genes + "</td>\n")
            sub_tree_display(path, results_url, subgroup_item.candidate_list, f, consider_homology, counter,
                             genes_names,
                             genes_list, data == 'removed_repetitions', genomeAssembly, use_thr)
    f.write(
        "<script>\nfunction hide(sgRNA) {\n\t var lst_r = document.getElementsByClassName(sgRNA.concat('_r'));\n    for(var i = 0; i < lst_r.length; ++i) {\n")
    f.write(
        "        if (lst_r[i].classList.contains('hideFromFilter')){lst_r[i].classList.remove('hideFromFilter');}\n        lst_r[i].style.display = 'none';\n    }\n    var lst = document.getElementsByClassName(sgRNA);\n    for(var i = 0; i < lst.length; ++i) {\n        lst[i].style.display = lst[i].style.display === 'none' ? 'table-cell' : 'none';")
    f.write("    }\n  }\n ")
    f.write(
        "function show_row(sgRNA) {\n\t var lst_r = document.getElementsByClassName(sgRNA.concat('_r'));\n \t var mm = document.getElementById(\"selectMM\").value;    for(var i = 0; i < lst_r.length; ++i) {\n")
    f.write(
        "        var myMM = lst_r[i].children[\"mms\"].innerText;\n \tif (parseInt(myMM)  <= parseInt(mm)) {\n lst_r[i].style.display = 'table-row';\n\t}\n \telse {lst_r[i].classList.add('hideFromFilter');}\n   }\n    var lst = document.getElementsByClassName(sgRNA);\n    for(var i = 0; i < lst.length; ++i) {\n        lst[i].style.display = lst[i].style.display === 'none' ? 'table-cell' : 'none';")
    f.write("    }\n  }\n ")
    f.write(
        "function show(elmnt,counter) {\n tmp = elmnt.innerText; \n if (tmp.startsWith(\"+\")) \n { \n elmnt.innerText = \"-\"+tmp.substring(1); \n} \n else {elmnt.innerText = \"+\"+tmp.substring(1);} \n document.getElementById(counter).style.display = document.getElementById(counter).style.display === 'none' ? 'block' : 'none'; \n updateQueryString(); \n}")
    f.write(
        "function filterMM(){\n var mm = document.getElementById(\"selectMM\").value;\n var tables = document.getElementsByClassName(\"table\");\n for (var t = 0; t < tables.length; ++t)\n {\n \ttable = tables[t];\n \ttr = table.getElementsByTagName(\"tr\");\n \tfor (i = 0; i < tr.length; i++)\n \t{\n    if(tr[i].style.display === 'none'){\n      if (tr[i].classList.contains('hideFromFilter')){}\nelse{continue;}\n    } \t\tvar row = tr[i].children[\"mms\"];\n \t\t if (row){\n \t\t\t var myMM = row.innerText;\n \t\t\t if (parseInt(myMM)  > parseInt(mm)) {\n        tr[i].classList.add('hideFromFilter');\n        tr[i].style.display = \"none\";\n       }\n \t\t\t else{tr[i].style.display = \"table-row\";tr[i].classList.remove('hideFromFilter');}\n \t\t}\n \t}\n }\n updateQueryString();\n}\n")

    if results_url:
        if not consider_homology:
            f.write(
                "function updateQueryString() { \n var mm = document.getElementById(\"selectMM\").value; \n var n = $(\"#scrl\").val(); \n document.getElementById(\"abc\").href=\"" + results_url + "&filterMM=\"+mm+\"&scrl=\"+n.toString();  \n }\n")
        else:
            f.write(
                "function updateQueryString() { \n var mm = document.getElementById(\"selectMM\").value; \n var innerSample= \"\"; \n innerSample += document.getElementById(1).style.display === 'block' ? \"1\" : \"0\"; \n innerSample += document.getElementById(2).style.display === 'block' ? \"1\" : \"0\"; \n innerSample += document.getElementById(3).style.display === 'block' ? \"1\" : \"0\"; \n innerSample += document.getElementById(4).style.display === 'block' ? \"1\" : \"0\"; \n var n = $(\"#scrl\").val(); \n document.getElementById(\"abc\").href=\"" + results_url + "&filterMM=\"+mm+\"&scrl=\"+n.toString()+\"&innerSample=\"+innerSample;  \n }\n")
        f.write(
            "$(window).scroll(function() { \n var n = $(window).scrollTop(); \n document.getElementById(\"scrl\").value = n.toString(); \n updateQueryString(); \n });\n")
    f.write("</script>\n")
    f.close()
