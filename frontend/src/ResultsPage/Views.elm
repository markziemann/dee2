module ResultsPage.Views exposing (..)

import Array
import Bool.Extra as BE
import Maybe.Extra as ME
import Dict
import Html exposing (..)
import Html.Attributes exposing (..)
import Html.Events as Events exposing (onClick)
import ResultsPage.Helpers exposing (..)
import ResultsPage.Types exposing (..)
import SearchPage.Types exposing (SearchResult, SearchResults)
import Set
import SharedTypes exposing (PaginationOffset)
import Table exposing (HtmlDetails, Status)


selectClickedResult : (SearchResult -> Msg) -> SelectedResults -> SearchResult -> List (Html.Attribute Msg)
selectClickedResult msg selectedResults searchResult =
    [ class (highlightRowIfTrue "table-primary" <| Dict.member searchResult.id selectedResults)
    , style "cursor" "pointer" -- First place not using bootstrap for style?
    , Events.onClick (msg searchResult)
    ]


stageResultForRemoval : (Int -> Msg) -> ResultsPendingRemoval -> SelectedResult -> List (Html.Attribute Msg)
stageResultForRemoval msg resultsPendingRemoval ( id, ( _, _ ) ) =
    [ class (highlightRowIfTrue "table-danger" <| Set.member id resultsPendingRemoval)
    , style "cursor" "pointer" -- First place not using bootstrap for style?
    , Events.onClick (msg id)
    ]


resultsTable : SelectedResults -> Table.Config SearchResult Msg
resultsTable selectedResults =
    Table.customConfig
        { toId = getId
        , toMsg = SetResultsTableState
        , columns =
            [ Table.stringColumn "Row" getId
            , Table.stringColumn "Species" (get "species")
            , Table.stringColumn "SRA Run" (get "SRR_accession")
            , Table.stringColumn "QC summary" (get "QC_summary")
            , Table.stringColumn "SRA experiment" (get "SRX_accession")
            , Table.stringColumn "SRA sample" (get "SRS_accession")
            , Table.stringColumn "SRA project" (get "SRP_accession")
            , Table.stringColumn "Sample" (get "Sample_name")
            , Table.stringColumn "Experiment" (get "GEO_series")
            ]
        , customizations =
            { defaultTable
                | tableAttrs =
                    [ class "table table-hover table-sm table-bordered table-responsive-lg"
                    , style "table-layout" "fixed" -- Prevents table going wider than parent element
                    , style "font-size" "clamp(12px, 4vw, 14px)"
                    ]
                , rowAttrs = selectClickedResult ResultClicked selectedResults
            }
        }


selectedTable : ResultsPendingRemoval -> Table.Config SelectedResult Msg
selectedTable resultsPendingRemoval =
    Table.customConfig
        { toId = Tuple.first >> String.fromInt
        , toMsg = SetSelectedResultsTableState
        , columns =
            [ Table.intColumn "Row" Tuple.first
            , Table.stringColumn "Species" (Tuple.second >> Tuple.first)
            , Table.stringColumn "SRA" (Tuple.second >> Tuple.second)
            ]
        , customizations =
            { defaultTable
                | tableAttrs =
                    [ class "table table-hover table-sm table-bordered"
                    , style "table-layout" "fixed" -- Prevents table going wider than parent element
                    , style "font-size" "clamp(12px, 4vw, 14px)"
                    ]
                , rowAttrs = stageResultForRemoval SelectedResultClicked resultsPendingRemoval
            }
        }


buttonOrSpinner : Bool -> SelectedResults -> Html Msg
buttonOrSpinner downloading rows =
    if not downloading then
        a
            [ hideWhenTrue "btn btn-outline-primary btn-block" (Dict.isEmpty rows) |> class
            , href <| "/download/" ++ queryString rows
            , attribute "download" "data.zip"
            , onClick DownloadRequested
            ]
            [ text "download" ]

    else
        button
            [ class "btn btn-primary btn-block"
            , attribute "type" "button"
            , attribute "disabled" "disabled"
            ]
            [ span
                [ class "spinner-border spinner-border-sm mr-2"
                , attribute "role" "status"
                , attribute "aria-hidden" "true"
                ]
                []
            , text "Loading..."
            ]


pageSelector disable page label =
    li [ class (disableWhenTrue "page-item" disable) ]
        [ button [ onClick (PageRequest page), class "page-link", attribute "tabindex" (BE.ifElse "-1" "0" disable) ]
            [ text label ]
        ]


pagination : SharedTypes.PaginationOffset -> Maybe ResultRows -> Html Msg
pagination ({ perPage, offset } as paginationOffset) maybeResultRows =
    let
        previousButton =
            if offset < perPage then
                pageSelector
                    True
                    paginationOffset
                    "Previous"

            else
                pageSelector False (PaginationOffset perPage <| offset - perPage) "Previous"

        resultRows = ME.unwrap Array.empty identity maybeResultRows

        nextButton =
            if Array.length resultRows < perPage then
                pageSelector
                    True
                    paginationOffset
                    "Next"

            else
                pageSelector
                    False
                    (PaginationOffset perPage <| offset + perPage)
                    "Next"
    in
    nav [ attribute "aria-label" "Page navigation" ]
        [ ul [ class "pagination justify-content-end" ]
            [ previousButton
            , nextButton
            ]
        ]


viewSearchResults : Model -> List (Html Msg)
viewSearchResults ({ searchResultRows, resultsTableState, resultsTableQuery, searchHits } as model) =
    [ div [ class "row" ]
        [ div [ class "col-xl-9" ]
            [ div [ class "bg-light text-primary" ]
                [ text "Hits: ", text (Maybe.withDefault "" (Maybe.map String.fromInt searchHits)) ]
            , Table.view
                (resultsTable model.selectedResults)
                resultsTableState
                (Maybe.withDefault [] (Maybe.map Array.toList searchResultRows))
            , pagination model.paginationOffset model.searchResultRows
            ]
        , div [ class "col-xl-3" ]
            [ div [ class "sticky-top" ]
                [ div [ class "bg-light text-primary" ]
                    [ text <| "Selected: " ++ (Dict.size model.selectedResults |> String.fromInt) ]
                , Table.view
                    (selectedTable model.resultsPendingRemoval)
                    model.selectedResultsTableState
                    (Dict.toList model.selectedResults)
                , div []
                    [ buttonOrSpinner model.downloading model.selectedResults
                    , button
                        [ hideWhenTrue
                            "btn btn-outline-danger btn-block"
                            (Set.isEmpty model.resultsPendingRemoval)
                            |> class
                        , onClick RemoveStagedSelections
                        ]
                        [ text "Remove" ]
                    ]
                ]
            ]
        ]
    ]
